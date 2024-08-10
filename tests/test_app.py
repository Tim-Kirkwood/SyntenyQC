# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:45:42 2024

@author: u03132tk
"""
from SyntenyQC.app import make_dirname, mixed_slashes, check_file_path_errors, check_dir_path_errors,read_args, check_args, main_cli
from general_mocks import mock_listdir, mock_get_gbk_files, mock_makedirs
import pytest

# =============================================================================
# Setup fixtures specific to this module and used across multiple TestClasses
# =============================================================================
@pytest.fixture
def setup_ok_file_path(monkeypatch):
    def mock_check_file_path_no_errors(long_var_name : str, path : str, 
                                       suffixes : list) -> tuple:
        return (None, None)
    monkeypatch.setattr('SyntenyQC.app.check_file_path_errors', 
                        mock_check_file_path_no_errors)
@pytest.fixture
def setup_ok_dir_path_with_gbk(monkeypatch):
    def mock_check_dir_path_no_errors(long_var_name : str, path : str) -> tuple:
        return (None, None)
    monkeypatch.setattr('SyntenyQC.app.check_dir_path_errors', 
                        mock_check_dir_path_no_errors)
    monkeypatch.setattr('SyntenyQC.app.get_gbk_files', 
                        mock_get_gbk_files)

@pytest.fixture
def setup_ok_dir_path_no_gbk(monkeypatch):
    def mock_check_dir_path_no_errors(long_var_name : str, path : str) -> tuple:
        return (None, None)
    def mock_get_no_gbk_files(folder : str) -> list:
        return []
    monkeypatch.setattr('SyntenyQC.app.get_gbk_files', 
                        mock_get_no_gbk_files)
    monkeypatch.setattr('SyntenyQC.app.check_dir_path_errors', 
                        mock_check_dir_path_no_errors)

@pytest.fixture 
def normal_setup(monkeypatch):
    def mock_isdir(path : str) -> bool:
        return True
    monkeypatch.setattr('os.path.isdir', 
                        mock_isdir)
    monkeypatch.setattr('os.listdir', 
                        mock_listdir) 
    #->"accession.gbk", "organism.gbk", "a_directory(1)", 'file3.gb', 'file3.txt', 'a_directory'
    monkeypatch.setattr('os.makedirs', 
                        mock_makedirs)
    #pass - dont make any dirs
    




# =============================================================================
# TestClasses
# =============================================================================
class TestMixedSlashes:
    
    @pytest.mark.parametrize('path,expected_bool', 
                             [('path/dir', False),
                              ('path\\dir', False),
                              ('path/to/dir', False),
                              ('path\\to\\dir', False),
                              ('path/to\\dir', True),
                              ('path\\to/dir', True)
                              ]
                             )                          
    def test_normal(self, path : str, expected_bool : bool):
        assert mixed_slashes(path) == expected_bool
    
    def test_exception(self):
        with pytest.raises(ValueError):
            mixed_slashes('path_with_no_slashes')

class TestMakeDirName:
        
    @pytest.fixture 
    def exception_setup(self, monkeypatch):
        def mock_is_not_dir(path : str) -> bool:
            return False
        monkeypatch.setattr('os.path.isdir', 
                            mock_is_not_dir)

        
    @pytest.mark.parametrize('base_dir_name,expected_name', 
                             [("accession", "a\\dir\\path\\accession"),
                              ("organism.gbk", "a\\dir\\path\\organism.gbk(1)"),
                              ("a_directory", "a\\dir\\path\\a_directory(2)"),
                              ("not_in_dir", "a\\dir\\path\\not_in_dir")
                              ]
                             )
    def test_normal(self, base_dir_name : str, expected_name : str, normal_setup):
        parent_dir = 'a\\dir\\path'
        assert make_dirname(parent_dir, base_dir_name) == expected_name
    
    def test_exception(self, exception_setup):
        with pytest.raises(ValueError):
            #raise valueerror if parent_dir does not exist
            make_dirname(parent_dir = 'not\\a\\parent\\dir', 
                         base_dir_name = 'base_dir_name')

class TestCheckFilePathErrors:
    
    @pytest.fixture
    def setup_isfile(self, monkeypatch):
        def mock_is_file(path : str) -> bool:
            return True
        monkeypatch.setattr('os.path.isfile', 
                            mock_is_file)
    
    @pytest.fixture
    def setup_isnotfile(self, monkeypatch): 
        def mock_is_not_file(path : str) -> bool:
            return False
        monkeypatch.setattr('os.path.isfile', 
                            mock_is_not_file)
    
    def test_normal(self, setup_isfile):
        code, message = check_file_path_errors('long_var_name', 
                                               path = 'path/with/slases.suffix', 
                                               suffixes = ['.suffix'])
        assert code == None
        assert message == None
    
    def test_path_fail(self, setup_isfile):
        code, message = check_file_path_errors('long_var_name', 
                                               path = 'path/with/mixed\\slashes.suffix', 
                                               suffixes = ['.suffix'])
        assert code == 1
        assert message == '--long_var_name path cannot contain forward and backward slashes.'
    
    @pytest.mark.parametrize('path, suffixes, expected_message',
                             [
                              ('path/with/slases.suffix',
                               ['.other_suffix'], 
                                "--long_var_name must be in csv format (as a file with suffix ['.other_suffix'])."),
                              ('path/with/slashes.suffix', 
                                ['suffix'], 
                                "--long_var_name must be in csv format (as a file with suffix ['suffix']).")
                              ]
                             )
    def test_suffix_fail(self, path : str, suffixes : list, expected_message : str, 
                         setup_isfile):
        code, message = check_file_path_errors('long_var_name', 
                                               path, 
                                               suffixes)
        assert code == 1
        assert message == expected_message
        
    
    def test_no_file(self, setup_isnotfile):
        code, message = check_file_path_errors('long_var_name', 
                                               path = 'path/with/slases.suffix', 
                                               suffixes = ['.suffix'])
        assert code == 2
        assert message == '--long_var_name file does not exist.'

class TestCheckDirPathErrors:
    @pytest.fixture
    def setup_isdir(self, monkeypatch):
        def mock_is_dir(path : str) -> bool:
            return True
        monkeypatch.setattr('os.path.isdir', 
                            mock_is_dir)
    
    @pytest.fixture
    def setup_isnotdir(self, monkeypatch):
        def mock_is_not_dir(path : str) -> bool:
            return False
        monkeypatch.setattr('os.path.isdir', 
                            mock_is_not_dir)
        
    def test_normal(self, setup_isdir):
        code, message = check_dir_path_errors('long_var_name', 
                                              path = 'path/with/slases')
        assert code == None
        assert message == None
    
    def test_path_fail(self, setup_isdir):
        code, message = check_dir_path_errors('long_var_name', 
                                              path = 'path/with/mixed\\slashes')
        assert code == 1
        assert message == '--long_var_name path cannot contain forward and backward slashes.'
    
    def test_no_dir(self, setup_isnotdir):
        code, message = check_dir_path_errors('long_var_name', 
                                              path = 'path/with/slases')
        assert code == 2
        assert message == '--long_var_name dir does not exist.'
    



class TestReadArgsCollect:
    
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
    def test_normal(self, command : str, setup_ok_file_path):
        #note - bp suffix ignored due to mock
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
                              'collect -gf a/folder -sf 0.5',
                              'collect -fake_flag fake_param',
                              'collect -fake_flag',
                              'collect -ns 1 -em an_email@domain.com',
                              'collect -ns fake_param', 
                              'collect -bp a/path -ns 1 -em an_email@domain.com -fn not_recognised -sp -wg'
                              ]
                             )
    def test_argparse_fail(self, cmd : str):
        with pytest.raises(SystemExit):
            args, parser = read_args(cmd.split())
        #TODO - check specific error message.  Need to clarify how argparse returns 
        #errors - capsys fixture returns captured_error with '\x08' special character 
        #at varios locations, making comparison difficult. 
                                   
    def test_no_binary_file(self, monkeypatch, capsys): 
        def mock_check_file_path_errors(long_var_name : str, path : str, 
                                        suffixes : list) -> tuple:
            return (2, 
                    f'--{long_var_name} file does not exist.')
        monkeypatch.setattr('SyntenyQC.app.check_file_path_errors',
                            mock_check_file_path_errors)
        cmd = 'collect -bp not/a/path.csv -ns 1 -em an_email@domain.com'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--binary_path file does not exist.'
        
    def test_suffix_error(self, monkeypatch, capsys):
        def mock_check_file_path_errors(long_var_name : str, path : str, 
                                        suffixes : list) -> tuple:
            return (1, 
                    f'--{long_var_name} must be in csv format (as a file with suffix {suffixes}).')
        monkeypatch.setattr('SyntenyQC.app.check_file_path_errors',
                            mock_check_file_path_errors)
        cmd = 'collect -bp a/path.not_csv -ns 1 -em an_email@domain.com'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == "--binary_path must be in csv format (as a file with suffix ['.txt', '.csv'])."
        
    def test_path_error(self, capsys):
        cmd = 'collect -bp a/mixed\\path.csv -ns 1 -em an_email@domain.com'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--binary_path path cannot contain forward and backward slashes.'
        
    def test_ns_error(self, setup_ok_file_path, capsys):
        cmd = 'collect -bp a/path -ns 0 -em an_email@domain.com'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--neighbourhood_size must be greater than 0'
    
    def test_email_error(self, setup_ok_file_path, capsys):
        cmd = 'collect -bp a/path -ns 1 -em an_email.domain.com'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--email is not correct - no @'


class TestReadArgsSieve:
    
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
    def test_normal(self, command : str, setup_ok_dir_path_with_gbk):
        args, parser = read_args(command.split())
        check_args(args, parser)
        assert args.command == 'sieve'
        assert args.genbank_folder == 'a\\folder'
        assert args.similarity_filter == 0.5
        
        if '-ev' in command:
            assert args.e_value == 0.05
        else:
            #default
            assert args.e_value == 10**-5
            
        if '-mi' in command:
            assert args.min_percent_identity == 60
        else:
            #default
            assert args.min_percent_identity == 50
        
        if '-mts' in command:
            assert args.max_target_seqs == 250
        else:
            #default
            assert args.max_target_seqs == 200
        
        if '-mev' in command:
            assert args.min_edge_view == 0.1
        else:
            #default
            assert args.min_edge_view == args.similarity_filter
            
        
    @pytest.mark.parametrize('cmd', 
                             ['sieve',
                              'sieve -bp a/path -ns 1 -em an_email@domain.com',
                              'sieve -fake_flag fake_param',
                              'sieve -sf 0.5',
                              'sieve -gf a/folder',
                              'sieve -gf fake_param', 
                              ]
                             )
    def test_argparse_fail(self, cmd : str):
        with pytest.raises(SystemExit):
            args, parser = read_args(cmd.split())    
        
        
    def test_no_dir(self, monkeypatch, capsys):
        def mock_check_dir_path_errors(long_var_name : str, path : str) -> tuple:
            return (2,
                    f'--{long_var_name} dir does not exist.')
        monkeypatch.setattr('SyntenyQC.app.check_dir_path_errors', 
                            mock_check_dir_path_errors)
        cmd = 'sieve -gf not/a/dir -sf 0.5'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--genbank_folder dir does not exist.'
    
    def test_path_error(self, capsys):
        cmd = 'sieve -gf a/mixed\\path -sf 0.5'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--genbank_folder path cannot contain forward and backward slashes.'
        
    def test_empty_dir(self, setup_ok_dir_path_no_gbk, capsys):
        cmd = 'sieve -gf empty/dir -sf 0.5'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == 'No genbank (.gbk, .gb) files in empty\\dir'
    
    @pytest.mark.parametrize('cmd', 
                             ['sieve -gf a/path -sf 0',
                              'sieve -gf a/path -sf 1.1'
                              ]
                             )
    def test_sf_error(self, cmd : str, setup_ok_dir_path_with_gbk, capsys):
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--similarity_filter must be between >0 and <=1.'
        
    @pytest.mark.parametrize('cmd', 
                             ['sieve -gf a/path -sf 0.1 -mev 0',
                              'sieve -gf a/path -sf 0.1 -mev 1.2',
                              ]
                             )
    def test_mev_error_range(self, cmd : str, setup_ok_dir_path_with_gbk, capsys):
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--min_edge_view must be between >0 and <=1.'
            
    def test_mev_error_greater_than_sf(self, setup_ok_dir_path_with_gbk, capsys):
        cmd = 'sieve -gf a/path -sf 0.1 -mev 0.2'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--min_edge_view must be <= similarity_filter.'
    
    @pytest.mark.parametrize('cmd', 
                             ['sieve -gf a/path -sf 0.5 -ev 0',
                              'sieve -gf a/path -sf 0.5 -ev 1.1',
                              ]
                             )
    def test_ev_error(self, cmd : str, setup_ok_dir_path_with_gbk, capsys):
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--e_value must be between >0 and <=1.'
    
    @pytest.mark.parametrize('cmd', 
                             ['sieve -gf a/path -sf 0.5 -mi 0',
                              'sieve -gf a/path -sf 0.5 -mi 101'
                              ]
                             )
    def test_mi_error(self, cmd : str, setup_ok_dir_path_with_gbk, capsys):
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--min_percent_identity must be between >0 and <=100.'
    
    def test_mts_error(self, setup_ok_dir_path_with_gbk, capsys):
        cmd = 'sieve -gf a/path -sf 0.5 -mts -1'
        args, parser = read_args(cmd.split())
        with pytest.raises(SystemExit):
            check_args(args, parser)
        captured = capsys.readouterr()
        assert captured.err == '--max_target_seqs must be between >0.' 
    
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
                             ['sieve', 'a\\folder', 10**-5, 50, 200, 0.5, 0.5, 
                              'a\\folder\\sieve_results']
                             ),
                           
                            ('sieve -gf a/folder -ev 0.05 -sf 0.5',
                             ['sieve', 'a\\folder', 0.05, 50, 200, 0.5, 0.5, 
                              'a\\folder\\sieve_results']
                             ),
                            
                            ('sieve -gf a/folder -ev 0.05 -mi 60 -sf 0.5',
                             ['sieve', 'a\\folder', 0.05, 60, 200, 0.5, 0.5, 
                              'a\\folder\\sieve_results']
                             ),
                           
                            ('sieve -gf a/folder -ev 0.05 -mi 60 -mts 250 -sf 0.5',
                             ['sieve', 'a\\folder', 0.05, 60, 250, 0.5, 0.5, 
                              'a\\folder\\sieve_results']
                             ),
                           
                            ('sieve -gf a/folder -ev 0.05 -mi 60 -mts 250 -mev 0.1 -sf 0.5',
                             ['sieve', 'a\\folder', 0.05, 60, 250, 0.5, 0.1, 
                              'a\\folder\\sieve_results']
                             )
                         ]
                        )
def test_main_cli(cmd : str, expected_params : list, monkeypatch, 
                  setup_ok_file_path, #for collect binary
                  setup_ok_dir_path_with_gbk, #for collect genbank folder
                  normal_setup, #for general file handling - isdir, listdir
                  capsys):

    # =============================================================================
    #     setup test specific monkeypatches
    # =============================================================================
    def mock_collect(binary_path : str, strict_span : bool, 
                     neighbourhood_size : int, write_genomes : bool, email : str, 
                     filenames : str, results_dir : str) -> str:
        nonlocal cli_params
        assert cli_params == [] #should only be set once per run
        cli_params = ['collect', binary_path, strict_span, neighbourhood_size, 
                      write_genomes, email, filenames, results_dir]
        return results_dir
    
    def mock_sieve(input_genbank_dir : str, e_value : float, 
                   min_percent_identity : int, max_target_seqs : int,
                   similarity_filter : float, results_dir : str, 
                   min_edge_view : float) -> str:
        nonlocal cli_params
        assert cli_params == []
        cli_params = ['sieve', input_genbank_dir, e_value, min_percent_identity, 
                      max_target_seqs, similarity_filter, min_edge_view, results_dir]
        return results_dir
    
    monkeypatch.setattr('SyntenyQC.app.collect', mock_collect)
    monkeypatch.setattr('SyntenyQC.app.sieve', mock_sieve)
    
    # =============================================================================
    #     run test
    # =============================================================================
    cli_params = []
    with pytest.raises(SystemExit):
        main_cli(cmd.split())
    assert cli_params == expected_params
    captured = capsys.readouterr()
    assert captured.err == f'Succsfully ran SyntenyQC -{cli_params[0]}.  Results in {cli_params[-1]}.'