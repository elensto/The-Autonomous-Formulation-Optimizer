"""Autonomous formulation optimizer: BO -> Experiment Agent/PhysCheck -> Simulation -> Result parser/Post-PhysCheck -> Observation store."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Make the team's original packages importable without modifying them.
for p in (ROOT / 'packages/BO/src', ROOT / 'packages/BO/reference', ROOT / 'packages/Agents/agents'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
