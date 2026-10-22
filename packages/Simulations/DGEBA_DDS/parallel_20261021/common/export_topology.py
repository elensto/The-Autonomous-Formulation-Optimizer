"""Remove only bond/react bookkeeping from a write_data export, preserving physics.

This is a fresh reactive stage, not a bitwise restart: stabilization/cooldown tags
are reset. Binary checkpoints remain available for exact stage continuation.
"""
from pathlib import Path
import hashlib
import json
import sys

def export(source, destination):
    content=source.read_bytes()
    lines=content.splitlines(keepends=True)
    allowed={'bond_react_props_internal','bond_react_stabilization_internal'}
    physical=[];removed={};section=None
    for line in lines:
        text=line.decode('utf-8').split('#')[0].strip()
        if text and (text[0].isalpha() or text[0]=='_'):
            name=text.split()[0]
            if name.startswith('bond_react_') and name not in allowed:
                raise ValueError('Unknown bond/react section: '+name)
            section=name if name in allowed else None
            if section:removed[section]=0
        if section:
            if text and text[0].isdigit():removed[section]+=1
        else:physical.append(line)
    if not removed:raise ValueError('Expected bond/react bookkeeping sections missing')
    result=b''.join(physical)
    destination.write_bytes(result)
    audit={'source':str(source),'destination':str(destination),
           'source_sha256':hashlib.sha256(content).hexdigest(),
           'destination_sha256':hashlib.sha256(result).hexdigest(),
           'removed_bookkeeping_rows':removed,
           'note':'Physical sections copied byte-for-byte. New stage resets reaction stabilization/cooldown tags.'}
    destination.with_suffix(destination.suffix+'.export.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2))

if __name__=='__main__':export(Path(sys.argv[1]),Path(sys.argv[2]))
