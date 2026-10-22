"""Create execution copies from verified originals; preserve all numerical data."""
from pathlib import Path
import hashlib
import json
import re

ROOT=Path(__file__).resolve().parents[1]
def main():
    raw=ROOT/'raw/lammps'
    source=raw/'data/small/data.100dgeba0_50dds_1gcm3_1'
    text=source.read_text()
    # The published small file lacks the separator after Dihedrals. LAMMPS
    # consumes the next line as a separator, silently dropping dihedral ID 1.
    adapted,n=re.subn(r'(?m)^Dihedrals[^\S\n]*\n(?=\S)','Dihedrals\n\n',text)
    assert text.split()==adapted.split(), 'Numerical/topological content changed'
    (ROOT/'inputs/data.initial').write_text(adapted,encoding='utf-8',newline='\n')
    for name in ['preTop.txt','postTop.txt','2preTop.txt','2postTop.txt','map1.txt','map2.txt']:
        s=(raw/'fixBondReact'/name).read_text().replace('BondingIDs','InitiatorIDs')
        (ROOT/'inputs'/name).write_text(s,encoding='utf-8',newline='\n')
    s=(raw/'inputScripts/minimizationEquilibration/in.min').read_text()
    a=s.index('mass       12');b=s.index('####\n',a)
    (ROOT/'inputs/coefficients_reaction.inc').write_text('# Source coefficients retained without numeric modifications.\n'+s[a:b],encoding='utf-8',newline='\n')
    changes={'source_file':str(source.relative_to(ROOT)),
             'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
             'execution_sha256':hashlib.sha256((ROOT/'inputs/data.initial').read_bytes()).hexdigest(),
             'blank_separator_inserted_after_Dihedrals':n==1,
             'source_and_execution_nonwhitespace_tokens_identical':True,
             'map_change':'BondingIDs -> InitiatorIDs; atom mappings retained',
             'coefficients_change':'Extracted reaction coefficients without numeric edits'}
    (ROOT/'reports/input_adaptations.json').write_text(json.dumps(changes,indent=2),encoding='utf-8')
    print(json.dumps(changes,indent=2))

if __name__=='__main__':main()
