import math
import unittest

from chinese_remainder_theorem import crt, CRTInputError


class TestCRT(unittest.TestCase):

    def test_single_congruence(self):
        # A single congruence is the trivial CRT case: the residue reduced mod the modulus.
        self.assertEqual(crt([2], [5]), 2)
        self.assertEqual(crt([7], [5]), 2)

    def test_classic_textbook_example(self):
        # x == 2 (mod 3), x == 3 (mod 5), x == 2 (mod 7)  =>  x == 23 (mod 105)
        self.assertEqual(crt([2, 3, 2], [3, 5, 7]), 23)

    def test_two_congruences(self):
        # x == 1 (mod 4), x == 2 (mod 9)  =>  x == 29 (mod 36)
        self.assertEqual(crt([1, 2], [4, 9]), 29)

    def test_result_in_range_and_correct(self):
        residues = [1, 4, 9]
        moduli = [5, 7, 11]
        product = math.prod(moduli)
        x = crt(residues, moduli)
        self.assertGreaterEqual(x, 0)
        self.assertLess(x, product)
        for r, m in zip(residues, moduli):
            self.assertEqual(x % m, r % m)

    def test_large_coprime_moduli(self):
        # Distinct primes are always pairwise coprime.
        moduli = [97, 101, 103, 107]
        residues = [10, 20, 30, 40]
        x = crt(residues, moduli)
        product = math.prod(moduli)
        self.assertLess(x, product)
        for r, m in zip(residues, moduli):
            self.assertEqual(x % m, r % m)

    def test_negative_residues_normalised(self):
        # x == -1 (mod 5) is x == 4 (mod 5); combine with x == 2 (mod 7).
        # x == 4 (mod 5), x == 2 (mod 7) => 4, 9, 14, 19, 24, 29 => 9 (mod 35)
        self.assertEqual(crt([-1, 2], [5, 7]), 9)

    def test_residue_larger_than_modulus_normalised(self):
        # x == 8 (mod 5) is x == 3 (mod 5); combine with x == 2 (mod 7).
        # x == 3 (mod 5), x == 2 (mod 7) => 3, 8, 13, 18, 23 (mod 35)
        self.assertEqual(crt([8, 2], [5, 7]), 23)

    def test_zero_modulus_rejected(self):
        with self.assertRaises(CRTInputError):
            crt([1, 2], [3, 0])

    def test_negative_modulus_rejected(self):
        with self.assertRaises(CRTInputError):
            crt([1, 2], [3, -5])

    def test_non_coprime_moduli_rejected(self):
        # gcd(6, 9) == 3, so the system may have no consistent solution.
        with self.assertRaises(CRTInputError):
            crt([1, 2], [6, 9])

    def test_mismatched_lengths_rejected(self):
        with self.assertRaises(CRTInputError):
            crt([1, 2, 3], [5, 7])
        with self.assertRaises(CRTInputError):
            crt([1], [5, 7])

    def test_empty_inputs_rejected(self):
        with self.assertRaises(CRTInputError):
            crt([], [])

    def test_non_integer_modulus_rejected(self):
        with self.assertRaises(CRTInputError):
            crt([1, 2], [3, 5.0])

    def test_non_integer_residue_rejected(self):
        with self.assertRaises(CRTInputError):
            crt([1, 2.0], [3, 5])

    def test_bool_modulus_rejected(self):
        # bool is a subclass of int; we reject it explicitly because modulus=True
        # is almost certainly a caller mistake rather than a genuine modulus of 1.
        with self.assertRaises(CRTInputError):
            crt([0], [True])

    def test_crt_input_error_is_value_error(self):
        self.assertTrue(issubclass(CRTInputError, ValueError))


if __name__ == "__main__":
    unittest.main()
