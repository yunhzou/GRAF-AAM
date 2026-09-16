"""Check that publication artifacts contain no execution-site information.

CPU model names, algorithm parameters, timings, and author affiliations remain
publishable. The check does not read or modify raw experiment archives.
"""
from pathlib import Path
import io
import json
import re

from pypdf import PdfReader

# Match concrete infrastructure identifiers, not chemistry partitions or CPUs.
PRIVATE_TEXT = re.compile(
    r"/(?:lustre|project|home|h|Users|scratch|mnt)/[^\s\"']+"
    r"|\b(?:cpu|gpu)-\d{4}\b"
    r"|\boci-[A-Za-z0-9_.-]+"
    r"|\b[A-Za-z0-9_.-]*(?:login-\d+)[A-Za-z0-9_.-]*\b"
    r"|\b(?:squeue|sacct|sbatch|slurm|nrt_enter|cpu-short)\b"
    r"|\b(?:concurrency limit|logical CPUs|GiB RAM|CPU nodes)\b",
    re.IGNORECASE,
)
PRIVATE_FIELDS = {
    'hostname', 'node_name', 'nodelist', 'slurm_job_id', 'slurm_partition',
    'cpus_per_task', 'mem_per_cpu', 'execution_host', 'login_host',
    'scheduler', 'allocation_id', 'cluster_name', 'concurrency',
}


def check_publication(name, payload):
    """Fail before packaging a PDF, document or evidence file with private fields."""
    suffix = Path(name).suffix.lower()
    if suffix == '.pdf':
        text = '\n'.join(p.extract_text() or '' for p in PdfReader(io.BytesIO(payload)).pages)
    elif suffix in {'.tex', '.json', '.md', '.svg', '.bib', '.html'}:
        text = payload.decode('utf-8')
    elif suffix in {'.py', '.sh'}:
        text = payload.decode('utf-8')
        # Code may name scheduler concepts, but must not embed site identifiers.
        concrete = re.compile(r'/(?:lustre|project|home|h|Users|scratch|mnt)/[^\s\"\']+|\b(?:cpu|gpu)-\d{4}\b|\boci-[A-Za-z0-9_.-]+', re.IGNORECASE)
        if concrete.search(text):
            raise ValueError(f'Execution-site identifier found in publication source: {name}')
        return
    elif suffix == '.sbatch':
        raise ValueError(f'Execution script is not a publication artifact: {name}')
    else:
        return
    if PRIVATE_TEXT.search(text):
        # Do not echo the sensitive matched value to logs.
        raise ValueError(f'Execution-site information found in publication artifact: {name}')
    if suffix == '.json':
        def walk(value):
            if isinstance(value, dict):
                if PRIVATE_FIELDS.intersection(k.lower() for k in value):
                    raise ValueError(f'Execution metadata field found in publication artifact: {name}')
                for child in value.values():
                    walk(child)
            elif isinstance(value, list):
                for child in value:
                    walk(child)
        walk(json.loads(text))
