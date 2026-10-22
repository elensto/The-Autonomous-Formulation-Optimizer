"""Fetch public sources, resuming interrupted Zenodo files and checking MD5."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
SOURCES=ROOT/'sources'
HEADERS={'User-Agent':'FormulationLab scientific reproducibility download'}

def open_url(url,extra=None):
    return urllib.request.urlopen(urllib.request.Request(url,headers=HEADERS | (extra or {})),timeout=90)

def main():
    SOURCES.mkdir(parents=True,exist_ok=True)
    (ROOT/'raw').mkdir(exist_ok=True)
    meta=SOURCES/'zenodo_record_7273800.json'
    if not meta.exists():
        with open_url('https://zenodo.org/api/records/7273800') as r:
            data=json.load(r)
        meta.write_text(json.dumps(data,indent=2),encoding='utf-8')
    data=json.loads(meta.read_text())
    for item in data['files']:
        dest=SOURCES/item['key']
        def valid():
            return dest.exists() and dest.stat().st_size==item['size'] and 'md5:'+hashlib.md5(dest.read_bytes()).hexdigest()==item['checksum']
        for attempt in range(5):
            if valid():break
            offset=dest.stat().st_size if dest.exists() else 0
            if offset>=item['size']:
                # A full-length invalid file cannot be repaired by appending.
                dest.rename(dest.with_name(dest.name+'.invalid-'+str(attempt)))
                offset=0
            try:
                with open_url(item['links']['self'],{'Range':f'bytes={offset}-'} if offset else None) as response:
                    append=offset>0 and response.status==206 and response.headers.get('Content-Range','').startswith(f'bytes {offset}-')
                    with dest.open('ab' if append else 'wb') as output:
                        shutil.copyfileobj(response,output)
            except (OSError,urllib.error.URLError) as error:
                print('Transfer interrupted; will resume:',dest.name,type(error).__name__,flush=True)
        if not valid():raise RuntimeError('Download incomplete or checksum mismatch: '+dest.name)
        print('Verified download',dest.name,flush=True)
    paper=SOURCES/'Livraghi_2023_Block_Chemistry.pdf'
    if not paper.exists():
        with open_url('https://fulir.irb.hr/8699/1/livraghi-et-al-2023-block-chemistry-for-accurate-modeling-of-epoxy-resins.pdf') as response:
            with paper.open('wb') as output:shutil.copyfileobj(response,output)
    assert paper.read_bytes().startswith(b'%PDF') and b'%%EOF' in paper.read_bytes()[-2048:]
    subprocess.run([sys.executable,str(ROOT/'scripts/verify_sources.py')],check=True)

if __name__=='__main__':main()
