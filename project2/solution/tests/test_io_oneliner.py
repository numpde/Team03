# RA, 2020-11-05

import io

from unittest import TestCase
from tempfile import NamedTemporaryFile


class TestOneliner(TestCase):
    def test_import(self):
        from idiva.io import Oneliner
        Oneliner

    def test_call_and_iterate(self):
        from idiva.io import Oneliner
        with NamedTemporaryFile(mode='w') as tf:
            with open(tf.name, mode='r') as fd:
                Oneliner(fd)

    def test_reads_lines(self):
        from idiva.io import Oneliner
        reference = ["Line1", "Line2"]
        with NamedTemporaryFile(mode='w') as tf:
            print(*reference, file=tf, sep='\n', flush=True)
            with open(tf.name, mode='r') as fd:
                candidate = list(Oneliner(fd))
                self.assertListEqual(reference, candidate)

    def test_skips_empty(self):
        from idiva.io import Oneliner
        reference = ["Line1", "Line2"]
        content = [""] + reference[:1] + [""] + reference[1:] + [""]
        with NamedTemporaryFile(mode='w') as tf:
            print(*content, file=tf, sep='\n', flush=True)
            with open(tf.name, mode='r') as fd:
                candidate = list(Oneliner(fd))
                self.assertListEqual(reference, candidate)

    def test_uninitialized_behavior(self):
        from idiva.io import Oneliner
        with NamedTemporaryFile(mode='w') as tf:
            with open(tf.name, mode='r') as fd:
                with self.assertRaises(AssertionError):
                    Oneliner(fd).last

    def test_keeps_last(self):
        from idiva.io import Oneliner
        reference = ["Line1", "Line2"]
        with NamedTemporaryFile(mode='w') as tf:
            print(*reference, file=tf, sep='\n', flush=True)
            with open(tf.name, mode='r') as fd:
                candidate = Oneliner(fd)
                self.assertEqual(reference[0], next(iter(candidate)))
                self.assertEqual(reference[0], candidate.last)
                self.assertEqual(reference[0], candidate.last)
                self.assertEqual(reference[1], next(iter(candidate)))
                self.assertEqual(reference[1], candidate.last)

    def test_reads_textio(self):
        from idiva.io import Oneliner
        reference = ["Line1", "Line2"]
        with NamedTemporaryFile(mode='w') as tf:
            print(*reference, file=tf, sep='\n', flush=True)
            with open(tf.name, mode='rb') as fd:
                with io.TextIOWrapper(fd) as fd:
                    candidate = list(Oneliner(fd))
                    self.assertListEqual(reference, candidate)

    def test_fails_binary(self):
        from idiva.io import Oneliner
        reference = ["Line1", "Line2"]
        with NamedTemporaryFile(mode='w') as tf:
            print(*reference, file=tf, sep='\n', flush=True)
            with open(tf.name, mode='rb') as fd:
                with self.assertRaises(TypeError):
                    list(Oneliner(fd))

