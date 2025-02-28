"""
RATTACA assignment: A tool for assigning HS West rats to RATTACA projects

This package provides functionality to assign rats to projects based on
genetic predictions, with the goal of maximizing between-group differences
for specific traits of interest.
"""

__version__ = '0.1.0'

from rattaca_assign.models.request import Request
from rattaca_assign.models.request_rattaca import RATTACA
from rattaca_assign.models.request_random import RandProj
from rattaca_assign.models.request_breeders import HSWBreeders

__all__ = ['Request', 'RATTACA', 'RandProj', 'HSWBreeders']