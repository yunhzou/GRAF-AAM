from setuptools import Extension, setup
import pybind11


setup(
    ext_modules=[
        Extension(
            "graft._group_ops",
            sources=["src/graft/native/group_ops.cpp"],
            include_dirs=[pybind11.get_include()],
            language="c++",
            extra_compile_args=["-O3", "-std=c++17"],
        ),
        Extension(
            "graft._native",
            sources=["src/graft/native/paired_mapping.cpp"],
            language="c++",
            extra_compile_args=["-O3", "-std=c++17"],
        )
    ]
)
