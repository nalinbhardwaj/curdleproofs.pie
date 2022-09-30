from typing import List
from py_ecc.secp256k1.secp256k1 import bytes_to_int
from py_ecc.optimized_bls12_381.optimized_curve import curve_order

from util import Fr

class CurdleproofsTranscript:
  def __init__(self) -> None:
    self.transcript: List[int] = []

  def append(self, label: bytes, item: bytes) -> None:
    self.transcript.append(label)
    self.transcript.append(item)
  
  def append_list(self, label: bytes, items: List[bytes]) -> None:
    self.transcript.append(label)
    for item in items:
      self.transcript.append(item)

  def get_and_append_challenge(self, label: bytes) -> Fr:
    # TODO: implement STROBE?
    self.transcript.append(label)
    buf = b''
    for i in range(0, 64):
      buf += bytes([i])
    return Fr(bytes_to_int(buf))
    # return Fr.one()

  def get_and_append_challenges(self, label: bytes, n: int) -> List[Fr]:
    return [self.get_and_append_challenge(label) for i in range(0, n)]

transcript = CurdleproofsTranscript()
transcript.append(b'hello', b'world')
# print(transcript.get_and_append_challenge(b'challenge'))
