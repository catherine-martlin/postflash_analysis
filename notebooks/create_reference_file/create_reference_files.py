import pandas as pd


def create_reference_file(year, working_directory, today, cadence, postflash_data):
    """This function will use input information to run the stack()
    function and create specific file names for the outputs.
    
    Parameters
    ----------
    year: int
        The year as an integer input.
    
    working_directory: str
        The path the files will be saved to. Needs trailing '\'
    
    today: str
        Format isn't important, needed to add date to filenames.
    
    cadence: int
        Number of years you want stacked. 1 = 1 year, 2 = biyearly, etc.
    
    postflash_data: pandas dataframe
        Dataframe to define the data you are using.
    
    Returns
    -------
    paths_year: str
        List of all FITS file needed to stack.
    
    outfile_year: str
        Filename for the outfile specified by year.
    
    error_outfile_year: str
        Filename for the error outfile specified by year.
    
    """
    fullframe_pf = postflash_data.loc[(postflash_data['subarray'] == False)]
    shutters = ('A', 'B')
    for shutter in shutters:
        fullframe_pf = postflash_data.loc[(postflash_data['subarray'] == False) & (postflash_data['shutter'] == '{}'.format(shutter)) & (postflash_data['flash_cur'] == 'MED') & (postflash_data['flash_dur'] == 100.0)]
        if cadence == 1:
            if year == 2012:
                fullframe_pf_year = fullframe_pf[(fullframe_pf['datetime'] > '{}-01-01 00:00:00'.format(str(year))) & (fullframe_pf['datetime'] < '{}-11-14 00:00:00'.format(str(year+1)))]
            else:
                fullframe_pf_year = fullframe_pf[(fullframe_pf['datetime'] > '{}-01-01 00:00:00'.format(str(year))) & (fullframe_pf['datetime'] < '{}-01-01 00:00:00'.format(str(year+1)))]
            paths_year = fullframe_pf_year.path.tolist()
            print(len(paths_year))
            outfile_year = '{}{}_fullframe_{}_flc_stack_{}.fits'.format(working_directory,str(year), shutter, today)
            error_outfile_year = '{}{}_fullframe_{}_flc_error_stack_{}.fits'.format(working_directory,str(year), shutter, today)
            print(outfile_year)
        else:
            fullframe_pf = postflash_data.loc[(postflash_data['subarray'] == False) & (postflash_data['shutter'] == '{}'.format(shutter)) & (postflash_data['flash_cur'] == 'MED') & (postflash_data['flash_dur'] == 100.0)]
            fullframe_pf_year = fullframe_pf[(fullframe_pf['datetime'] > '{}-01-01 00:00:00'.format(str(year))) & (fullframe_pf['datetime'] < '{}-01-01 00:00:00'.format(str(year+cadence)))]
            paths_year = fullframe_pf_year.path.tolist()
            print(len(paths_year))
            outfile_year = '{}{}_cadence{}_fullframe_{}_flc_stack_{}.fits'.format(working_directory,str(year), str(cadence), shutter, today)
            error_outfile_year = '{}{}_cadence{}_fullframe_{}_flc_error_stack_{}.fits'.format(working_directory,str(year), str(cadence), shutter, today)
            print(outfile_year)

    return paths_year, outfile_year, error_outfile_year, fullframe_pf_year
