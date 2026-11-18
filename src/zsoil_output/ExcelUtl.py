import xlrd
import xlwt

#====================================================
def EXCEL_get_sheets (xls_file):
# ====================================================
    wbk = xlrd.open_workbook (xls_file)
    return wbk.sheet_names()

#====================================================
def EXCEL_get_data_from_sheet (xls_file,sheet_index):
#====================================================
    wbk = xlrd.open_workbook (xls_file)
    sheet = wbk.sheet_by_index (sheet_index)
    #sheet = wbk.sheet_by_name (sheet_name)  #another format

    num_cols = sheet.ncols   # Number of columns
    data = []
    for row_index in range(0,sheet.nrows):    # Iterate through rows
        row_data = []
        for col_index in range(0,len(sheet.row(row_index))):  # Iterate through columns
            cell_obj = sheet.cell (row_index,col_index)
            row_data.append (cell_obj.value)
        data.append (row_data)
    return data


#====================================================
def EXCEL_get_data_from_sheet_by_name (xls_file,sheet_name):
#====================================================
    #print "OPEN WORKBOOK",sheet_name
    wbk = xlrd.open_workbook (xls_file)
    #sheet = wbk.sheet_by_index (sheet_index)
    sheet = wbk.sheet_by_name (sheet_name)  #another format

    num_cols = sheet.ncols   # Number of columns
    data = []
    for row_index in range(0,sheet.nrows):    # Iterate through rows
        row_data = []
        for col_index in range(0,len(sheet.row(row_index))):  # Iterate through columns
            cell_obj = sheet.cell (row_index,col_index)
            row_data.append (cell_obj.value)
        data.append (row_data)
    return data


def str_to_float (text):
    try:
        d = float(text)
    except ValueError:
        d = None
    return d

# ================================================
def ExtractColumnOfFloats (array,col):
# ================================================
    out = []
    for row in array:
        a = str_to_float ( row [col] )
        if a != None:
            out.append (float(a))
    return out


# ================================================
def ExtractColumnsOfFloats (array,cols):
# ================================================
    out = []
    for row in array:
        aux = []
        for k in range(len(cols)):
            a = str_to_float ( row [cols [k]] )
            if a != None:
                aux.append (float(a))
        if len(aux) == len(cols):
            out.append(aux)
    return out


def EXCEL_new_workbook ():
    wbk = xlwt.Workbook ()
    return wbk

def EXCEL_new_sheet (wbk,sheet_label):
    sheet = wbk.add_sheet (sheet_label)
    return sheet

def EXCEL_save_workbook (wbk,xls_file_name):
    wbk.save (xls_file_name)

def EXCEL_store_table (sheet,table,offs_rows=0,offs_cols=0):
    for i in range (len(table)):
        for j in range (len(table[i])):
            sheet.write (i+offs_rows,j+offs_cols,table [i][j])

def EXCEL_store_np_vector (sheet,vector,offs_rows=0,offs_cols=0):
    for i in range (vector.size):
    	sheet.write (i+offs_rows,offs_cols,vector [i])


def main ():
    tab = EXCEL_get_data_from_sheet ("1A.xls",0)
    depth = ExtractColumnOfFloats (tab,0)
    Mz    = ExtractColumnOfFloats (tab,4)
    print (depth)
    print (Mz)

if __name__ == '__main__':
    main()