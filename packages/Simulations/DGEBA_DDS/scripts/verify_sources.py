"""Verify downloaded Zenodo files before safe extraction and provenance recording."""
from pathlib import Path
import hashlib
import json
import zipfile

ROOT=Path(__file__).resolve().parents[1]
def main():
    meta=json.loads((ROOT/'sources/zenodo_record_7273800.json').read_text())
    records=[]
    for item in meta['files']:
        path=ROOT/'sources'/item['key']
        raw=path.read_bytes()
        assert len(raw)==item['size'], f'Truncated file: {path.name}'
        assert 'md5:'+hashlib.md5(raw).hexdigest()==item['checksum'], f'Checksum mismatch: {path.name}'
        records.append({'file':path.name,'bytes':len(raw),'source_url':item['links']['self'],
                        'zenodo_checksum':item['checksum'],'sha256':hashlib.sha256(raw).hexdigest()})
        if zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as archive:
                assert archive.testzip() is None
                for entry in archive.infolist():
                    assert (ROOT/'raw'/entry.filename).resolve().is_relative_to((ROOT/'raw').resolve())
                archive.extractall(ROOT/'raw')
        print('Verified',path.name)
    paper=ROOT/'sources/Livraghi_2023_Block_Chemistry.pdf'
    assert paper.read_bytes().startswith(b'%PDF')
    records.append({'file':paper.name,'sha256':hashlib.sha256(paper.read_bytes()).hexdigest(),
                    'source_url':'https://fulir.irb.hr/8699/1/livraghi-et-al-2023-block-chemistry-for-accurate-modeling-of-epoxy-resins.pdf'})
    (ROOT/'sources/download_manifest.json').write_text(json.dumps({'record':7273800,'doi':'10.1021/acs.jpcb.3c04724','files':records},indent=2),encoding='utf-8')

if __name__=='__main__':main()
