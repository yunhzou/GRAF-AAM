const fs=require('fs');const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
(async()=>{
 const browser=await chromium.launch({executablePath:process.env.CHROME_PATH,headless:true,args:['--enable-unsafe-swiftshader']});
 const page=await browser.newPage({viewport:{width:1600,height:1100}}),errors=[],requests=[];
 page.on('pageerror',e=>errors.push(String(e)));page.on('request',r=>{if(/^https?:/.test(r.url()))requests.push(r.url())});
 await page.goto('file://'+process.argv[2]);await page.waitForFunction(()=>typeof CASES!=='undefined');
 const count=await page.locator('#list .case').count();const checks=[];if(count!==140)throw Error('Expected all 140 cases');
 for(let i=0;i<count;i++){
  await page.locator('#list .case').nth(i).click();const frame=page.frameLocator('#viewer');
  await frame.locator('#interpPlayBtn').waitFor();await frame.locator('#interp_method').filter({hasText:'persistent_fragment'}).waitFor();
  const btns=frame.locator('#mech-sel button');const choices=await btns.count();
  for(let j=0;j<choices;j++){
   await btns.nth(j).click();await frame.locator('#interpFrame').fill('50');await frame.locator('#interpFrame').dispatchEvent('input');
   await frame.locator('#showClashHighlights').check();
   if(!(await frame.locator('#interp_t').textContent()).includes('0.50'))throw Error('Slider failed');
  }
  await frame.locator('#interpPlayBtn').click();await page.waitForTimeout(220);await frame.locator('#interpPlayBtn').click();
  if((await frame.locator('#interp_t').textContent()).includes('0.50'))throw Error('Playback failed');
  checks.push({casePosition:i,choices});if(i%20===0)console.log('Checked '+(i+1)+'/'+count);
 }
 await page.locator('#list .case').nth(3).click();await page.frameLocator('#viewer').locator('#interpPlayBtn').waitFor();await page.evaluate(()=>scrollTo(0,0));await page.screenshot({path:process.argv[3]+'.png'});if(errors.length||requests.length)throw Error(JSON.stringify({errors,requests}));
 fs.writeFileSync(process.argv[3]+'.json',JSON.stringify({count,checks,errors,requests},null,2));await browser.close();console.log(JSON.stringify({count,checks}));
})().catch(e=>{console.error(e);process.exit(1)});
