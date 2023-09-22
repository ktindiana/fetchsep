#!/usr/bin/env python

"""Tests for `fetchsep` package."""


import unittest

import fetchsep


class TestFetchsep(unittest.TestCase):
    """Tests for `fetchsep` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_version(self):
        """Test that version exists."""
        self.assertIsInstance(fetchsep.__version__, str)
