# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:45:42 2024

@author: u03132tk
"""
#import argparse
#import SyntenyQC.app
from SyntenyQC.app import make_dirname, mixed_slashes, check_file_path_errors, check_dir_path_errors,read_args, check_args, main_cli
#import os
from general_mocks import mock_listdir
import pytest

#import os.path

#print(os.path.abspath(SyntenyQC.app.__file__))
#raise ValueError
#https://stackoverflow.com/a/18161115/11357695

@pytest.mark.parametrize('base_dir_name,expected_name', 
                         [("accession", "a\\dir\\path\\accession"),
                          ("organism.gbk", "a\\dir\\path\\organism.gbk(1)"),
                          ("a_directory", "a\\dir\\path\\a_directory(2)"),
                          ("not_in_dir", "a\\dir\\path\\not_in_dir")
                          ]
                         )
def test_make_dirname(base_dir_name, expected_name, monkeypatch):
    def mock_isdir(path):
        return True
    monkeypatch.setattr('os.path.isdir', mock_isdir)
    monkeypatch.setattr('os.listdir', mock_listdir)
    assert make_dirname('a\\dir\\path', base_dir_name) == expected_name

def test_make_dirname_exception(monkeypatch):
    def mock_is_not_dir(path):
        return False
    monkeypatch.setattr('os.path.isdir', mock_is_not_dir)
    with pytest.raises(ValueError):
        make_dirname('a\\parent\\dir', 'base_dir_name')

@pytest.mark.parametrize('path,expected_bool', 
                         [('path/dir', False),
                          ('path\\dir', False),
                          ('path/to/dir', False),
                          ('path\\to\\dir', False),
                          ('path/to\\dir', True),
                          ('path\\to/dir', True)
                          ]
                         )                          
def test_mixed_slashes(path, expected_bool):
    assert mixed_slashes(path) == expected_bool

def test_mixed_slashes_fail():
    with pytest.raises(ValueError):
        mixed_slashes('path_with_no_slashes')

@pytest.mark.parametrize('path,mock_func, expected_code, expected_message',
                         [('path/with/mixed\\slashes.suffix',       
                          None, 
                          1, 
                          '--long_var_name path cannot contain forward and backward slashes.'),
                         ('path/with/slases.suffix',       
                          'is_dir', 
                          None,
                          None),
                         ('path/with/slases.suffix',       
                          'is_not_dir', 
                          2, 
                          '--long_var_name dir does not exist.')
                         ])
def test_check_dir_path_errors(path,mock_func, expected_code, expected_message, monkeypatch):
    def mock_is_dir(path):
        return True
    def mock_is_not_dir(path):
        return False
    if mock_func == 'is_dir': 
        monkeypatch.setattr('os.path.isdir', mock_is_dir)
    else:
        monkeypatch.setattr('os.path.isdir', mock_is_not_dir)

    
    code, message = check_dir_path_errors('long_var_name', 
                                      path)
    assert code == expected_code
    assert message == expected_message

    
@pytest.mark.parametrize('path,suffixes,mock_func, expected_code, expected_message',
                         [('path/with/mixed\\slashes.suffix', 
                           ['other_suffix'], 
                           None, 
                           1, 
                           '--long_var_name path cannot contain forward and backward slashes.'),
                          ('path/with/slases.suffix', 
                           ['.suffix'], 
                           'is_file', 
                           None,
                           None),
                          ('path/with/slases.suffix',
                           ['other_suffix'],
                           'is_file', 
                            1, 
                            "--long_var_name must be in csv format (as a file with suffix ['other_suffix'])."),
                          ('path/with/slases.suffix', 
                           ['other_suufix'], 
                           'is_not_file',
                           2,
                           '--long_var_name file does not exist.')
                          ]
                         )
def test_check_file_path_errors(path,suffixes,mock_func, expected_code, expected_message, monkeypatch):
    def mock_is_not_file(path):
        return False
    def mock_is_file(path):
        return True
    
    if mock_func == 'is_file': 
        monkeypatch.setattr('os.path.isfile', mock_is_file)
    else:
        monkeypatch.setattr('os.path.isfile', mock_is_not_file)
    code, message = check_file_path_errors('long_var_name', 
                                      path, 
                                      suffixes)
    assert code == expected_code
    assert message == expected_message

def mock_check_file_path_no_errors(long_var_name, path, suffixes):
    return (None, None)

def mock_check_dir_path_no_errors(long_var_name, path):
    return (None, None)

def mock_check_dir_path_errors(long_var_name, path):
    return (2,f'--{long_var_name} dir does not exist.')
def mock_check_file_path_errors(long_var_name, path, suffixes):
    return (2, f'--{long_var_name} file does not exist.')

@pytest.mark.parametrize('command', 
                         ['collect -bp a/path -ns 1 -em an_email@domain.com',     
                          'collect -bp a/path -ns 1 -em an_email@domain.com -fn organism',
                          'collect -bp a/path -ns 1 -em an_email@domain.com -fn organism -sp',
                          'collect -bp a/path -ns 1 -em an_email@domain.com -fn organism -sp -wg',
                          'collect -bp a/path -ns 1 -em an_email@domain.com -fn accession',
                          'collect -bp a/path -ns 1 -em an_email@domain.com -fn accession -sp',
                          'collect -bp a/path -ns 1 -em an_email@domain.com -fn accession -sp -wg',
                          'collect -bp a\\path -ns 1 -em an_email@domain.com',     
                          'collect -bp a\\path -ns 1 -em an_email@domain.com -fn organism',
                          'collect -bp a\\path -ns 1 -em an_email@domain.com -fn organism -sp',
                          'collect -bp a\\path -ns 1 -em an_email@domain.com -fn organism -sp -wg',
                          'collect -bp a\\path -ns 1 -em an_email@domain.com -fn accession',
                          'collect -bp a\\path -ns 1 -em an_email@domain.com -fn accession -sp',
                          'collect -bp a\\path -ns 1 -em an_email@domain.com -fn accession -sp -wg'
                          ]
                         )
def test_arg_parsing_collect(command, monkeypatch):
    monkeypatch.setattr('SyntenyQC.app.check_file_path_errors', 
                        mock_check_file_path_no_errors)
    args, parser = read_args(command.split())
    check_args(args, parser)
    assert args.command == 'collect'
    assert args.binary_path == 'a\\path'
    assert args.neighbourhood_size == 1
    assert args.email == 'an_email@domain.com'
    if '-fn accession' in command:
        assert args.filenames == 'accession'
    else:
        #default is organism
        assert args.filenames == 'organism'
    if '-sp' in command:
        assert args.strict_span == True
    else:
        assert args.strict_span == False
    if '-wg' in command:
        assert args.write_genomes == True
    else:
        assert args.write_genomes == False

@pytest.mark.parametrize('cmd', 
                         ['collect',
                          'collect -fake_flag fake_param',
                          'collect -fake_flag',
                          'collect -ns 1 -em an_email@domain.com',
                          'collect -ns fake_param', 
                          'collect -bp a/path -ns 1 -em an_email@domain.com -fn not_recognised -sp -wg', 
                          
                          ]
                         )
def test_arg_parsing_collect_read_fail(cmd, monkeypatch):
    with pytest.raises(SystemExit):
        args, parser = read_args(cmd.split())
    #TODO - clarify how argparse returns errors 
    #the following approach isn't working atm - messages broadly match but for some reason, captured_error gets '\x08'
    #special character at varios locations that makes comparison difficult. Also applies for -h.
    # captured = capsys.readouterr()
    
    # if cmd == 'collect -bp a/path -ns 1 -em an_email@domain.com -fn not_recognised -sp -wg':
    #     specific_error = "SyntenyQC collect: error: argument -fn/--filenames: "\
    #                       "invalid choice: 'not_recognised' (choose from "\
    #                       "'organism', 'accession')\n"
    # elif cmd == 'collect -ns fake_param':
    #     specific_error = "SyntenyQC collect: error: argument -ns/--neighbourhood_size: invalid int value: 'fake_param'\n"
    # else:
    #     specific_error = 'SyntenyQC collect: error: the following arguments are '\
    #                       'required: -bp/--binary_path, -ns/--neighbourhood_size, '\
    #                       '-em/--email\n'
    # expected_message = 'usage: SyntenyQC collect [-h] -bp \x08 -ns \x08 -em \x08 [-fn ] [-sp] [-wg]\n' + specific_error
    
    # print (f'CAPTURED ERROR:\n[{captured.err}]')
    # print (expected_message)
    # print (captured.err == expected_message) 
    # print ([i if i == ii else (i,ii) for i, ii in zip(expected_message, captured.err) ])    
    # assert captured.err == expected_message
                               


                                       
@pytest.mark.parametrize('cmd,message',                                        
                         [('collect -bp not/a/path -ns 1 -em an_email@domain.com', 
                           '--binary_path file does not exist.'),     
                          ('collect -bp a/path -ns -1 -em an_email@domain.com -fn organism',
                           '--neighbourhood_size must be greater than 0'),
                          ('collect -bp a/mixed\\path -ns 1 -em an_email@domain.com -fn organism -sp',
                           '--binary_path path cannot contain forward and backward slashes.'),
                          ('collect -bp a/path -ns 1 -em an_email.domain.com -fn organism -sp -wg',
                           '--email is not correct - no @')
                          ]
                         )#TODO check more path types
def test_arg_parsing_collect_check_fail(cmd, message, monkeypatch, capsys):
    if '-bp not/a/path' in cmd:
        monkeypatch.setattr('SyntenyQC.app.check_file_path_errors', mock_check_file_path_errors)
    elif '-bp a/mixed\\path' in cmd:
        pass
    else:
        monkeypatch.setattr('SyntenyQC.app.check_file_path_errors', mock_check_file_path_no_errors)
        
    args, parser = read_args(cmd.split())
    with pytest.raises(SystemExit):
        check_args(args, parser)
    captured = capsys.readouterr()
    print (captured.err)
    assert captured.err == message
    
    
    
    
    
@pytest.mark.parametrize('command', 
                         
                         ['sieve -gf a/folder -sf 0.5',
                         'sieve -gf a/folder -ev 0.05 -sf 0.5',
                         'sieve -gf a/folder -ev 0.05 -mi 60 -sf 0.5',
                         'sieve -gf a/folder -ev 0.05 -mi 60 -mts 250 -sf 0.5',
                         'sieve -gf a/folder -ev 0.05 -mi 60 -mts 250 -mev 0.1 -sf 0.5',
                         'sieve -gf a\\folder -sf 0.5',
                         'sieve -gf a\\folder -ev 0.05 -sf 0.5',
                         'sieve -gf a\\folder -ev 0.05 -mi 60 -sf 0.5',
                         'sieve -gf a\\folder -ev 0.05 -mi 60 -mts 250 -sf 0.5',
                         'sieve -gf a\\folder -ev 0.05 -mi 60 -mts 250 -mev 0.1 -sf 0.5',
                         ]
                         
                         )
def test_arg_parsing_sieve(command, monkeypatch):
    monkeypatch.setattr('SyntenyQC.app.check_dir_path_errors', mock_check_dir_path_no_errors)
    args, parser = read_args(command.split())
    check_args(args, parser)
    assert args.command == 'sieve'
    assert args.genbank_folder == 'a\\folder'
    assert args.similarity_filter == 0.5
    
    if '-ev' in command:
        assert args.e_value == 0.05
    else:
        assert args.e_value == 10**-5
        
    if '-mi' in command:
        assert args.min_percent_identity == 60
    else:
        assert args.min_percent_identity == 50
    
    if '-mts' in command:
        assert args.max_target_seqs == 250
    else:
        assert args.max_target_seqs == 200
    
    if '-mev' in command:
        assert args.min_edge_view == 0.1
    else:
        assert args.min_edge_view == args.similarity_filter
        
    
@pytest.mark.parametrize('cmd', 
                         ['sieve',
                          'sieve -fake_flag fake_param',
                          'sieve -sf 0.5',
                          'sieve -gf a/folder',
                          'sieve -gf fake_param', 
                          
                          ]
                         )
def test_arg_parsing_sieve_read_fail(cmd, monkeypatch):
    with pytest.raises(SystemExit):
        args, parser = read_args(cmd.split())    
    
    
@pytest.mark.parametrize('cmd,message',                                        
                         [('sieve -gf not/a/dir -sf 0.5', 
                           '--genbank_folder dir does not exist.'),     
                          ('sieve -gf a/mixed\\path -sf 0.5',
                           '--genbank_folder path cannot contain forward and backward slashes.'),
                          
                          ('sieve -gf a/path -sf 0', 
                           '--similarity_filter must be between >0 and <=1.'),
                          ('sieve -gf a/path -sf 1.1', 
                           '--similarity_filter must be between >0 and <=1.'),
                          
                          ('sieve -gf a/path -sf 0.1 -mev 0.2', 
                           '--min_edge_view must be <= similarity_filter.'),
                          ('sieve -gf a/path -sf 0.1 -mev -0.2', 
                           '--min_edge_view must be between >0 and <=1.'),
                          ('sieve -gf a/path -sf 0.1 -mev 1.2', 
                           '--min_edge_view must be between >0 and <=1.'),
                          
                          ('sieve -gf a/path -sf 0.5 -ev 0', 
                           '--e_value must be between >0 and <=1.'),
                          ('sieve -gf a/path -sf 0.5 -ev 1.1', 
                           '--e_value must be between >0 and <=1.'),
                          
                          ('sieve -gf a/path -sf 0.5 -mi 0', 
                           '--min_percent_identity must be between >0 and <=100.'),
                          ('sieve -gf a/path -sf 0.5 -mi 101', 
                           '--min_percent_identity must be between >0 and <=100.'),
                          
                          ('sieve -gf a/path -sf 0.5 -mts -1', 
                           '--max_target_seqs must be between >0.')
                          ]
                         )#TODO check more path types
def test_arg_parsing_sieve_check_fail(cmd, message, monkeypatch, capsys):
    if '-gf not/a/dir' in cmd:
        monkeypatch.setattr('SyntenyQC.app.check_dir_path_errors', mock_check_dir_path_errors)
    elif '-gf a/mixed\\path' in cmd:
        pass
    else:
        monkeypatch.setattr('SyntenyQC.app.check_dir_path_errors', mock_check_dir_path_no_errors)
        
    args, parser = read_args(cmd.split())
    with pytest.raises(SystemExit):
        check_args(args, parser)
    captured = capsys.readouterr()
    #print (captured.err)
    assert captured.err == message
    
    
@pytest.mark.parametrize('cmd,expected_params',
                         [
                             ('collect -bp a/path.suffix -ns 10 -em an_email@domain.com',  
                              ['collect', 'a\\path.suffix', False, 10, False, 
                               'an_email@domain.com', 'organism', 'a\\path']
                              ),
                                                                   
                             ('collect -bp a/path.suffix -ns 10 -em an_email@domain.com -fn organism',
                              ['collect', 'a\\path.suffix', False, 10, False, 
                               'an_email@domain.com', 'organism', 'a\\path']
                              ),
                          
                            ('collect -bp a/path.suffix -ns 10 -em an_email@domain.com -fn organism -sp',
                             ['collect', 'a\\path.suffix', True, 10, False, 
                              'an_email@domain.com', 'organism', 'a\\path']
                             ),
                           
                            ('collect -bp a/path.suffix -ns 10 -em an_email@domain.com -fn organism -sp -wg',
                             ['collect', 'a\\path.suffix', True, 10, True, 
                              'an_email@domain.com', 'organism', 'a\\path']
                             ),
                             
                            ('collect -bp a/path.suffix -ns 10 -em an_email@domain.com -fn accession',
                             ['collect', 'a\\path.suffix', False, 10, False, 
                              'an_email@domain.com', 'accession', 'a\\path']
                             ),
                            
                            ('collect -bp a/path.suffix -ns 10 -em an_email@domain.com -fn accession -sp',
                             ['collect', 'a\\path.suffix', True, 10, False, 
                              'an_email@domain.com', 'accession', 'a\\path']
                             ),
                           
                            ('collect -bp a/path.suffix -ns 10 -em an_email@domain.com -fn accession -sp -wg',
                             ['collect', 'a\\path.suffix', True, 10, True, 
                              'an_email@domain.com', 'accession', 'a\\path']
                             ),
                          

                              
                            ('sieve -gf a/folder -sf 0.5',
                             ['sieve', 'a\\folder', 10**-5, 50, 200, 0.5, 0.5, 'a\\folder\\ClusterSieve']
                             ),
                           
                            ('sieve -gf a/folder -ev 0.05 -sf 0.5',
                             ['sieve', 'a\\folder', 0.05, 50, 200, 0.5, 0.5, 'a\\folder\\ClusterSieve']
                             ),
                            
                            ('sieve -gf a/folder -ev 0.05 -mi 60 -sf 0.5',
                             ['sieve', 'a\\folder', 0.05, 60, 200, 0.5, 0.5, 'a\\folder\\ClusterSieve']
                             ),
                           
                            ('sieve -gf a/folder -ev 0.05 -mi 60 -mts 250 -sf 0.5',
                             ['sieve', 'a\\folder', 0.05, 60, 250, 0.5, 0.5, 'a\\folder\\ClusterSieve']
                             ),
                           
                            ('sieve -gf a/folder -ev 0.05 -mi 60 -mts 250 -mev 0.1 -sf 0.5',
                             ['sieve', 'a\\folder', 0.05, 60, 250, 0.5, 0.1, 'a\\folder\\ClusterSieve']
                             )
                         ]
                        )
def test_main_cli(cmd,expected_params,monkeypatch, capsys):
    def mock_makedirs(path):
        pass
    def mock_isdir(path):
        return True
    def mock_collect(binary_path, strict_span, neighbourhood_size, write_genomes, 
                     email, filenames, results_dir):
        nonlocal cli_params
        assert cli_params == None #should only be set once per run
        cli_params = ['collect', binary_path, strict_span, neighbourhood_size, 
                      write_genomes, email, filenames, results_dir]
        return results_dir
    def mock_sieve(genbank_folder, e_value, min_percent_identity, max_target_seqs,
                   similarity_filter, results_dir, min_edge_view):
        nonlocal cli_params
        assert cli_params == None
        cli_params = ['sieve', genbank_folder, e_value, min_percent_identity, 
                      max_target_seqs, similarity_filter, min_edge_view, results_dir]
        return results_dir
    monkeypatch.setattr('os.makedirs', mock_makedirs)
    monkeypatch.setattr('os.listdir', mock_listdir)
    monkeypatch.setattr('os.path.isdir', mock_isdir)
    monkeypatch.setattr('SyntenyQC.app.check_file_path_errors', 
                        mock_check_file_path_no_errors)
    monkeypatch.setattr('SyntenyQC.app.check_dir_path_errors', 
                        mock_check_dir_path_no_errors)
    monkeypatch.setattr('SyntenyQC.app.collect', mock_collect)
    monkeypatch.setattr('SyntenyQC.app.sieve', mock_sieve)
    cli_params = None
    with pytest.raises(SystemExit):
        main_cli(cmd.split())
    assert cli_params == expected_params
    captured = capsys.readouterr()
    assert captured.err == f'Succsfully ran SyntenyQC -{cli_params[0]}.  Results in {cli_params[-1]}.'