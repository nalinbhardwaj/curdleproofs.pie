from typing import List
from py_ecc.secp256k1.secp256k1 import bytes_to_int
from py_ecc.optimized_bls12_381.optimized_curve import curve_order
from merlin import MerlinTranscript

from util import Fr


class CurdleproofsTranscript(MerlinTranscript):
    def append(self, label: bytes, item: bytes) -> None:
        self.append_message(label, item)

    def append_list(self, label: bytes, items: List[bytes]) -> None:
        for item in items:
            self.append_message(label, item)

    def get_and_append_challenge(self, label: bytes) -> Fr:
        while True:
            challenge_bytes = self.challenge_bytes(label, 255)
            f = Fr(bytes_to_int(challenge_bytes))
            if f != Fr.zero():
                self.append(label, challenge_bytes)
                return f

    def get_and_append_challenges(self, label: bytes, n: int) -> List[Fr]:
        return [self.get_and_append_challenge(label) for _ in range(0, n)]


# transcript = CurdleproofsTranscript(b'curdleproofs')
# transcript.append(b'hello', b'world')
# print(transcript.get_and_append_challenge(b'challenge'))
