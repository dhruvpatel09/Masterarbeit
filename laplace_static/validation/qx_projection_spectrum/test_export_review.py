#!/usr/bin/env python3
"""Regression checks for large collection exports and partial-ZIP recovery."""
import tempfile
import unittest
import zipfile
from pathlib import Path
from export_review import export_archive, SMALL_LIMIT


class ExportReviewTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.small = self.root/'case.out'
        self.small.write_bytes(b'validated measurement')
        self.raw = self.root/'profiles_raw.csv'
        self.raw.write_bytes(b'x'*(SMALL_LIMIT+1))
        self.archive = self.root/'review.zip'

    def test_large_validated_collection_is_included(self):
        export_archive(self.archive, [self.small, self.raw], self.root, {self.raw})
        with zipfile.ZipFile(self.archive) as bundle:
            self.assertEqual(set(bundle.namelist()), {'case.out', 'profiles_raw.csv'})
            self.assertEqual(bundle.read('profiles_raw.csv'), self.raw.read_bytes())

    def test_unapproved_large_file_is_rejected_before_creation(self):
        with self.assertRaises(ValueError):
            export_archive(self.archive, [self.small, self.raw], self.root, set())
        self.assertFalse(self.archive.exists())

    def test_matching_partial_archive_can_be_completed_and_reused(self):
        with zipfile.ZipFile(self.archive, 'w') as bundle:
            bundle.write(self.small, 'case.out')
        self.assertEqual(export_archive(self.archive, [self.small, self.raw], self.root, {self.raw}), 'complete')
        original = self.archive.read_bytes()
        self.assertEqual(export_archive(self.archive, [self.small, self.raw], self.root, {self.raw}), 'already complete')
        self.assertEqual(original, self.archive.read_bytes())

    def test_changed_partial_archive_is_preserved(self):
        with zipfile.ZipFile(self.archive, 'w') as bundle:
            bundle.writestr('case.out', b'unrelated content')
        original = self.archive.read_bytes()
        with self.assertRaises(ValueError):
            export_archive(self.archive, [self.small, self.raw], self.root, {self.raw})
        self.assertEqual(original, self.archive.read_bytes())


if __name__ == '__main__':
    unittest.main()
