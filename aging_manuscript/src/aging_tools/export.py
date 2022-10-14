import datetime
import os
from aging_tools import inout

p_material = inout.get_internal_path('figures/tstoeger/material')


def get_material_path(p, insert_date_time=False):
    """
    Takes extension path p, and makes it as a subfile
    within material folder
    """

    p = os.path.join(
        p_material,
        p)

    if insert_date_time:
        [fo, fn] = os.path.split(p)
        [fb, ext] = os.path.splitext(fn)
        dt = datetime.datetime.today().strftime('%y%m%d_%H%M')
        p = os.path.join(fo, fb + '_' + dt + ext)

    # inout.ensure_presence_of_directory(p)

    return p


def export_image(p, insert_date_time=True):
    """
    Will export a current figure;
        - makes font edit-able in illustrator

    Input:
        insert_date_time    optional; default True: will insert
                                date, hour, minute before file
                                extension
    """

    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42  # edit-able in illustrator

    p = os.path.join(
        p_material,
        p)

    if insert_date_time:
        [fo, fn] = os.path.split(p)
        [fb, ext] = os.path.splitext(fn)
        dt = datetime.datetime.today().strftime('%y%m%d_%H%M')
        p = os.path.join(fo, fb + '_' + dt + ext)

    inout.ensure_presence_of_directory(p)

    mpl.pyplot.savefig(p, bbox_inches='tight')


def export_raster_image(p, dpi, insert_date_time=True):
    """
    Will export a current figure with 600 dpi

    Input:
        p       str path to file
        dpi     int dots per inch
        insert_date_time    optional; default True: will insert
                                date, hour, minute before file
                                extension
    """

    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42  # edit-able in illustrator

    p = get_material_path(p, insert_date_time)
    inout.ensure_presence_of_directory(p)

    mpl.pyplot.savefig(p, dpi=dpi, bbox_inches='tight')


def export_full_frame(p, df, insert_date_time=True, save_index=True):
    """
    Will export a dataframe to materials
        - makes font edit-able in illustrator

    Input:
        p                   subpath within material folder
        df                  dataframe to export
        insert_date_time    optional; default True: will insert
                                date, hour, minute before file
                                extension
        save_index          optional; default True: will also
                                export the index
    """

    p = os.path.join(
        p_material,
        p)

    if p.endswith('.csv.gz'):
        p = p[:-3]
        compress = True
        file_format = 'csv'
    elif p.endswith('.csv'):
        compress = False
        file_format = 'csv'
    elif p.endswith('.xlsx'):
        file_format = 'xlsx'
    else:
        raise EnvironmentError(
            'No support for preseent file type.')

    if insert_date_time:
        [fo, fn] = os.path.split(p)
        [fb, ext] = os.path.splitext(fn)
        dt = datetime.datetime.today().strftime('%y%m%d_%H%M')
        p = os.path.join(fo, fb + '_' + dt + ext)

    inout.ensure_presence_of_directory(p)

    if file_format == 'csv':
        if compress:
            p = p + '.gz'
            df.to_csv(p, compression='gzip', index=save_index)
        else:
            df.to_csv(p, index=save_index)
    elif file_format == 'xlsx':
        df.to_excel(p, index=save_index, engine='openpyxl')
