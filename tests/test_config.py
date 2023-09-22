#!/usr/bin/env python

"""Tests for `fetchsep.utils.config` package."""


import unittest

import fetchsep.utils.config as cfg


class TestConfig(unittest.TestCase):
    """Tests for `fetchsep.utils.config` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_datadir(self):
        """Test that version exists."""
        self.assertIsInstance(cfg.datapath, str)
