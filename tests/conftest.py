# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 11:32:25 2024

@author: u03132tk
"""
#from Bio import SeqIO
import logging
import pytest

@pytest.fixture 
def log_setup(monkeypatch):
    def mock_get_log_handlers(path : str) -> list:
        return [logging.StreamHandler()]
    monkeypatch.setattr('SyntenyQC.helpers.get_log_handlers', 
                        mock_get_log_handlers)
    


