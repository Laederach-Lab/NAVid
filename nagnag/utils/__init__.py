''' Classes for data/analysis on biological modules
* Genomic sequences (`Seq`)
* NAGNAG splice scenarios (`Scenario`)
* NAGNAG-affecting variants (`Variant`)
'''

__all__ = ["Seq", "Scenario", "Variant"]

from .sequence import Seq
from .scenarios import Scenario
from .variants import Variant
