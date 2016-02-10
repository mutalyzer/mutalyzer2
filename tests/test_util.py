"""
Tests for the mutalyzer.util module.
"""


from __future__ import unicode_literals

import pytest

from mutalyzer import util


@pytest.mark.parametrize('ref,var,descr,first,last_ref,last_var', [
    ('MTAPQQMT*', 'MTAQQMT*', 'p.(Pro4del)', 3, 4, 3),
    ('MTAPQQMT*', 'MTAQMT*', 'p.(Pro4_Gln5del)', 3, 5, 3),
    ('MTAPQQT*', 'MTAQQMT*', 'p.(Pro4_Gln6delinsGlnGlnMet)', 3, 6, 6),
    ('MTAPQQMT*', 'MTAPQQMTMQ*', 'p.(*9Metext*2)', 8, 9, 11),
    ('MTAPQQMT*', 'MTAPQQMTMQ', 'p.(*9Metext*?)', 8, 9, 10)])
def test_in_frame_description(ref, var, descr, first, last_ref, last_var):
    """
    In-frame description of difference between two proteins.
    """
    assert util.in_frame_description(ref, var) == (
        descr, first, last_ref, last_var)


@pytest.mark.parametrize('ref,var,descr,first,last_ref,last_var', [
    ('MTAPQQMT*', 'MTAQQMT*', 'p.(Pro4Glnfs*5)', 3, 9, 8),
    ('MTAPQQMT*', 'MTAQMT*', 'p.(Pro4Glnfs*4)', 3, 9, 7),
    ('MTAPQQT*', 'MTAQQMT*', 'p.(Pro4Glnfs*5)', 3, 8, 8),
    ('MTAPQQT*', 'MTAQQMT', 'p.(Pro4Glnfs*?)', 3, 8, 7)])
def test_out_of_frame_description(ref, var, descr, first, last_ref, last_var):
    """
    Out-of-frame description of difference between two proteins.
    """
    assert util.out_of_frame_description(ref, var) == (
        descr, first, last_ref, last_var)
