"""
Utility functions for writing excel sheets
"""


def safe_len(x):
    try:
        return len(x)
    except TypeError:
        return 16


def series_to_col_width(series):

    if series.dtype == 'float64':
        width = 16
    elif series.dtype == 'int64':
        width = 10
    else:
        width = int(series.map(safe_len).max()*1.3)

    return max(width, len(str(series.name)) + 3)


def write_sheet(df, sheet_name, writer):
    """
    Writes the sheet, but auto-formats the column width
    """
    df.to_excel(writer, sheet_name=sheet_name)
    worksheet = writer.sheets[sheet_name]
    df = df.reset_index()
    for idx in range(df.shape[1]):
        series = df.iloc[:, idx]
        max_len = series_to_col_width(series)
        worksheet.set_column(idx, idx, max_len)
