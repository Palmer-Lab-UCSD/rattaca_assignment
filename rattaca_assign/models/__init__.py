'''
Models (request class objects) for the rattaca_assign package.
'''

from rattaca_assign.models.request import Request
from rattaca_assign.models.request_rattaca import RATTACA
from rattaca_assign.models.request_random import RandProj
from rattaca_assign.models.request_breeders import HSWBreeders

__all__ = [
    'Request',
    'RATTACA',
    'RandProj',
    'HSWBreeders'
]