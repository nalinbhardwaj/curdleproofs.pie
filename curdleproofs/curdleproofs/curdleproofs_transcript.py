from typing import List
from merlin_transcripts import MerlinTranscript
from py_arkworks_bls12381 import Scalar
from curdleproofs.util import CURVE_ORDER


class CurdleproofsTranscript(MerlinTranscript):
    def append(self, label: bytes, item: bytes) -> None:
        self.append_message(label, item)

    def append_list(self, label: bytes, items: List[bytes]) -> None:
        for item in items:
            self.append_message(label, item)

    def get_and_append_challenge(self, label: bytes) -> Scalar:
        while True:
            challenge_bytes = self.challenge_bytes(label, 32)
            challenge_int = int.from_bytes(challenge_bytes, byteorder='little')
            if challenge_int >= CURVE_ORDER:
                continue
            # `Scalar.from_le_bytes` will error if value if >= CURVE_ORDER
            f = Scalar.from_le_bytes(challenge_bytes)
            if not f.is_zero():
                self.append(label, challenge_bytes)
                return f

    def get_and_append_challenges(self, label: bytes, n: int) -> List[Scalar]:
        return [self.get_and_append_challenge(label) for _ in range(0, n)]


# transcript = CurdleproofsTranscript(b'curdleproofs')
# transcript.append(b'hello', b'world')
# print(transcript.get_and_append_challenge(b'challenge'))
