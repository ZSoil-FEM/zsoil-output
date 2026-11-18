"""Read tabular project data from Microsoft Excel workbooks.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The module provides workbook, worksheet, cell-range, and conversion helpers used to import
auxiliary data for ZSoil post-processing workflows.
"""

import xlrd
import xlwt


def EXCEL_get_sheets(xls_file):
    """Open an Excel workbook and return its worksheets.

    :param xls_file: Xls file used by this operation.
    :type xls_file: str
    :return: Worksheet names in the workbook.
    :rtype: list[str]
    """
    wbk = xlrd.open_workbook(xls_file)
    return wbk.sheet_names()


def EXCEL_get_data_from_sheet(xls_file, sheet_index):
    """Return all cell values from a worksheet selected by index.

    :param xls_file: Xls file used by this operation.
    :type xls_file: str
    :param sheet_index: Sheet index used by this operation.
    :type sheet_index: int
    :return: All cell values from a worksheet selected by index.
    :rtype: list
    """
    wbk = xlrd.open_workbook(xls_file)
    sheet = wbk.sheet_by_index(sheet_index)
    # sheet = wbk.sheet_by_name (sheet_name)  #another format

    num_cols = sheet.ncols   # Number of columns
    data = []
    for row_index in range(0, sheet.nrows):    # Iterate through rows
        row_data = []
        for col_index in range(0, len(sheet.row(row_index))):  # Iterate through columns
            cell_obj = sheet.cell(row_index, col_index)
            row_data.append(cell_obj.value)
        data.append(row_data)
    return data


def EXCEL_get_data_from_sheet_by_name(xls_file, sheet_name):
    """Return all cell values from a worksheet selected by name.

    :param xls_file: Xls file used by this operation.
    :type xls_file: str
    :param sheet_name: Sheet name used by this operation.
    :type sheet_name: str
    :return: All cell values from a worksheet selected by name.
    :rtype: list
    """
    #print "OPEN WORKBOOK",sheet_name
    wbk = xlrd.open_workbook(xls_file)
    #sheet = wbk.sheet_by_index (sheet_index)
    sheet = wbk.sheet_by_name(sheet_name)  # another format

    num_cols = sheet.ncols   # Number of columns
    data = []
    for row_index in range(0, sheet.nrows):    # Iterate through rows
        row_data = []
        for col_index in range(0, len(sheet.row(row_index))):  # Iterate through columns
            cell_obj = sheet.cell(row_index, col_index)
            row_data.append(cell_obj.value)
        data.append(row_data)
    return data


def str_to_float(text):
    """Convert text to a floating-point value, returning zero for empty input.

    :param text: Text used by this operation.
    :type text: str
    :return: Converted floating-point value, or zero when the input is empty.
    :rtype: Optional[float]
    """
    try:
        d = float(text)
    except ValueError:
        d = None
    return d


def ExtractColumnOfFloats(array, col):
    """Return one table column converted to floating-point values.

    :param array: Array used by this operation.
    :type array: Sequence
    :param col: Column used by this operation.
    :type col: int
    :return: One table column converted to floating-point values.
    :rtype: list
    """
    out = []
    for row in array:
        a = str_to_float(row[col])
        if a != None:
            out.append(float(a))
    return out


def ExtractColumnsOfFloats(array, cols):
    """Return selected table columns converted to floating-point values.

    :param array: Array used by this operation.
    :type array: Sequence
    :param cols: Cols used by this operation.
    :type cols: Sequence[int]
    :return: The selected table columns converted to floating-point values.
    :rtype: list
    """
    out = []
    for row in array:
        aux = []
        for k in range(len(cols)):
            a = str_to_float(row[cols[k]])
            if a != None:
                aux.append(float(a))
        if len(aux) == len(cols):
            out.append(aux)
    return out


def EXCEL_new_workbook():
    """Create and return an empty Excel workbook.

    :return: New empty workbook.
    :rtype: xlwt.Workbook
    """
    wbk = xlwt.Workbook()
    return wbk


def EXCEL_new_sheet(wbk, sheet_label):
    """Add a named worksheet to an Excel workbook and return it.

    :param wbk: Wbk used by this operation.
    :type wbk: workbook
    :param sheet_label: Sheet label used by this operation.
    :type sheet_label: str
    :return: Newly added worksheet.
    :rtype: xlwt.Worksheet
    """
    sheet = wbk.add_sheet(sheet_label)
    return sheet


def EXCEL_save_workbook(wbk, xls_file_name):
    """Save an Excel workbook to the requested file.

    :param wbk: Wbk used by this operation.
    :type wbk: workbook
    :param xls_file_name: Xls file name used by this operation.
    :type xls_file_name: str
    """
    wbk.save(xls_file_name)


def EXCEL_store_table(sheet, table, offs_rows=0, offs_cols=0):
    """Write a two-dimensional table into a worksheet at the requested offset.

    :param sheet: Sheet used by this operation.
    :type sheet: worksheet
    :param table: Table used by this operation.
    :type table: Sequence[Sequence]
    :param offs_rows: Offs rows used by this operation.
    :type offs_rows: int
    :param offs_cols: Offs cols used by this operation.
    :type offs_cols: int
    """
    for i in range(len(table)):
        for j in range(len(table[i])):
            sheet.write(i+offs_rows, j+offs_cols, table[i][j])


def EXCEL_store_np_vector(sheet, vector, offs_rows=0, offs_cols=0):
    """Write a NumPy vector into a worksheet at the requested offset.

    :param sheet: Sheet used by this operation.
    :type sheet: worksheet
    :param vector: Vector used by this operation.
    :type vector: numpy.ndarray
    :param offs_rows: Offs rows used by this operation.
    :type offs_rows: int
    :param offs_cols: Offs cols used by this operation.
    :type offs_cols: int
    """
    for i in range(vector.size):
        sheet.write(i+offs_rows, offs_cols, vector[i])


def main():
    """Run this module's example or command-line workflow."""
    tab = EXCEL_get_data_from_sheet("1A.xls", 0)
    depth = ExtractColumnOfFloats(tab, 0)
    Mz = ExtractColumnOfFloats(tab, 4)
    print (depth)
    print (Mz)


if __name__ == '__main__':
    main()
