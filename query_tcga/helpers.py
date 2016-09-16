
from query_tcga.log_with import log_with

@log_with()
def _compute_start_given_page(page, size):
    """ compute start / from position given page & size
    """
    return (page*size+1)

@log_with()
def _convert_to_list(x):
    """ Convert x to a list if not already a list

    Examples
    -----------

    >>> _convert_to_list('Clinical')
    ['Clinical']
    >>> _convert_to_list(['Clinical'])
    ['Clinical']
    >>> _convert_to_list(('Clinical','Biospecimen'))
    ['Clinical', 'Biospecimen']

    """
    if not(x):
        return(None)
    elif isinstance(x, list):
        return(x)
    elif isinstance(x, str):
        return([x])
    else:
        return(list(x))
