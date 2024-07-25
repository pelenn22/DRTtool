import tkinter as tk
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import pathlib
import tksheet
import os


from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
from galvani import BioLogic



from scipy.optimize import nnls
from scipy.optimize import lsq_linear
from scipy.signal import find_peaks
from scipy.integrate import simps
from lmfit.models import GaussianModel
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from RangeSlider.RangeSlider import RangeSliderH, RangeSliderV
#from dataphile.statistics.regression.modeling import Parameter, Model, CompositeModel, AutoGUI
#from dataphile.statistics.distributions import gaussian1D
from datetime import datetime
from decimal import Decimal
from time import sleep


root = tk.Tk()
s = ttk.Style()
s.theme_use('winnative')
root.title("DRTtool")
root.iconbitmap("icon_drt.ico")

# bug fix for tkinter's treeview and colored lines
if root.getvar('tk_patchLevel') == '8.6.9':  # and OS_Name=='nt':
    def fixed_map(option):
        return [elm for elm in s.map('Treeview', query_opt=option) if elm[:2] != ('!disabled', '!selected')]
    s.map('Treeview', foreground=fixed_map('foreground'), background=fixed_map('background'))

# root.iconbitmap('tests/icon_drt.ico')
# root.iconphoto(True, tk.PhotoImage(file=r'C:\Users\peter\PycharmProjects\DataImport_test\tests\icon3n.png'))
#root.geometry("+0+50")
root.state('zoomed')
ttkstyle = ttk.Style(root)
ttkstyle.configure('lefttab.TNotebook', tabposition='nw')
my_notebook = ttk.Notebook(root, style='lefttab.TNotebook')

importtab = Frame(my_notebook, width=200, height=200, bg="grey")
optitab = Frame(my_notebook, bg="light grey")
fittab = Frame(my_notebook, bg="light grey")


my_notebook.add(optitab, text="DRT")
my_notebook.add(fittab, text="Fitting tool")
my_notebook.pack(expand=True, fill="both")

my_notebook.select(optitab)


rlabel = Label(root)
clabel = Label(root)
llabel = Label(root)
text_list = []
datav_list_all = []
origin_list = []
newdf = []
freqq = []
negimq = []
realsq = []
canvas = FigureCanvasTkAgg(Figure())
toolbar = NavigationToolbar2Tk(canvas, Frame())
fig_frame = Frame()
freq = []
negim = []
reals = []
df = []
df_list =[]
rtimes =[]
matCreg = []
bvecreg = []
n_result =[]
data1_list = []
xdata_g = []
ydata_g = []
entryset2 = []
descv = []
areav = []
thickv = []
tempv = []
symcheckvarv = 1
canvas3 = FigureCanvasTkAgg(Figure(), master=fittab)
canvas5 = FigureCanvasTkAgg(Figure(), master=fittab)
instlist = []
instlist2 = []
labellist_all = []
collist = ["Column 1",
           "Column 2",
           "Column 3",
           "Column 4",
           "Column 5",
           "Column 6",
           "Column 7"]
clickede1 = StringVar()
clickede1.set(collist[0])
clickede2 = StringVar()
clickede2.set(collist[1])
clickede3 = StringVar()
clickede3.set(collist[2])

headerv = IntVar()
headerv.set(12)
header0 = "None"
endv = IntVar()
endv.set(140)


freqcol = 0
realcol = 4
imagcol = 5

delimlist = [",", ";", "\s+", "\t"]
delimvalue = StringVar()
delimvalue.set(delimlist[0])

presetlist = ["Zahner", "Autolab", "Biologic", "Circuit Simulator"]
presetvalue = StringVar()



entryset = [0.001, 10, 30]
lamalist = np.geomspace(entryset[0], entryset[1], entryset[2])

arean = 1
lowerblog = -6
middleblog = -2

def multiimport2():
    global text_list
    global datav_list_all
    global origin_list
    global df_list
    global instlist
    global labellist_all
    global delimiterv
    if not import_tree2.get_children():
        pass
    else:
        import_tree2.delete(*import_tree2.get_children())
    text_files = filedialog.askopenfilenames()
    #print(text_files)
    text_list = list(text_files)

    count = 0
    labellist_all = []
    for record in text_list:
        path = pathlib.PurePath(record)
        import_tree2.insert(parent='', index='end', iid=count, text=path.name, values=record)
        labellist_all.append(path.name[:-4])
        count += 1

    origin_list = []
    x = import_tree2.get_children()
    child_id = import_tree2.get_children()[0]
    import_tree2.focus(child_id)
    import_tree2.selection_set(child_id)


    #item_text = [import_tree2.item(item_i, "values") for item_i in x]
    item_text = [text_files[i] for i in range(len(text_files))]
    #print(item_text)
    #print(item_text[0][0][-4:])
    if item_text[0][0][-4:] == '.mpr':
        files_all = [BioLogic.MPRfile(file) for file in item_text]
        df_all = [pd.DataFrame(files_all[i].data) for i in range(len(files_all))]
        datav_list_all = [data_list_i.values.tolist()[0:] for data_list_i in df_all]
    else:
        delimv0 = delimvalue.get()
        headv = headerv.get()
        data_list = [pd.read_csv(file, sep=delimv0, header=headv, encoding='latin1') for file in item_text]
        datav_list_all = [data_list_i.values.tolist()[0:] for data_list_i in data_list]

    #import_settings()
    instlist = []
    create_instlists()
    import_tree2.bind('<ButtonRelease-1>', plotselection2)
    plotselection2(0)

def multiimport_drt():
    global text_list2
    global datav_list_all2
    global origin_list2
    global instlist2
    global labellist_all2
    global delimiterv

    if not import_drt_tree.get_children():
        pass
    else:
        import_drt_tree.delete(*import_drt_tree.get_children())
    text_files2 = filedialog.askopenfilenames()
    text_list2 = list(text_files2)

    count = 0
    labellist_all2 = []
    for record in text_list2:
        path = pathlib.PurePath(record)
        import_drt_tree.insert(parent='', index='end', iid=count, text=path.name, values=record)
        labellist_all2.append(path.name[:-4])
        count += 1

    origin_list2 = []
    x = import_drt_tree.get_children()
    child_id = import_drt_tree.get_children()[0]
    import_drt_tree.focus(child_id)
    import_drt_tree.selection_set(child_id)

    item_text = [import_drt_tree.item(item_i, "values") for item_i in x]
    data_list = [pd.read_csv(file[0], sep='\t', header=0, encoding='latin1') for file in item_text]
    datav_list_all2 = [data_list_i.values.tolist() for data_list_i in data_list]

    instlist2 = []
    create_instlists2()
    import_drt_tree.bind('<ButtonRelease-1>', plotselection3)
    plotselection3(0)
    #plotselection2(0)

def create_canvas2(masterframe, size):
    global canvas2
    global toolbar
    global ax1
    global ax2
    global ax3
    global ax2t
    global figure
    (w, h) = size
    inchsize = (w / 25.4, h / 25.4)
    figure = Figure(inchsize)
    ax1 = figure.add_subplot(211, aspect='equal', xlabel='Re(Z)/\u03a9 cm²',
                             ylabel='-Im(Z)/\u03a9 cm²')  # Nyquist Plot
    ax2 = figure.add_subplot(212, xscale='log', xlabel='log(tau/s)',
                             ylabel='res')  # drt plot
    ax2t = ax2.secondary_xaxis('top')
    ax2t.set_xlabel('f/Hz')
    ax2t.set_xscale('log')
    ax2t.xaxis.set_major_formatter(FuncFormatter(lambda x, _: '$10^{{{}}}$'.format(int(np.log10(1 / 2*np.pi*x)))))
    ax2.set_ylabel("h(tau)/Ohm cm2")
    #ax3.set_ylabel("imaginary residual /%")
    # ax3.set_ylabel("imag residues")

    # create canvas as matplotlib drawing area
    canvas2 = FigureCanvasTkAgg(figure, master=masterframe)
    toolbar = NavigationToolbar2Tk(canvas2, masterframe)
    canvas2.get_tk_widget().pack_forget()
    canvas2.get_tk_widget().pack()
    figure.tight_layout()

def create_canvas4(masterframe, size):
    global canvas4
    global toolbar4
    global figurek
    global ax4
    global ax5

    (w, h) = size
    inchsize = (w / 25.4, h / 25.4)
    figurek = Figure(inchsize)
    # self.ax1 = self.figure.add_subplot(211, aspect='equal', xlabel='Re(Z)/Omega',
    # ylabel='-Im(Z)/Omega')  # Nyquist Plot
    ax4 = figurek.add_subplot(111, xscale='log', xlabel='Frequency/Hz',
                                       ylabel='res')  # reals reisdual plot
    ax5 = ax4.twinx()  # imaginary residual plot
    ax4.set_ylabel("real residues")
    ax5.set_ylabel("imag residues", color='C1')
    ax5.spines['right'].set_color('C1')
    ax5.tick_params(axis='y', colors='C1', which='both')
    ax4.format_coord = lambda x, y: ''
    ax5.format_coord = lambda x, y: ''
    # create canvas as matplotlib drawing area
    canvas4 = FigureCanvasTkAgg(figurek, master=masterframe)
    toolbar4 = NavigationToolbar2Tk(canvas4, masterframe)
    canvas4.get_tk_widget().pack_forget()
    canvas4.get_tk_widget().pack()
    figurek.tight_layout()

def create_canvas_fit(masterframe, size):
    global canvas_fit
    global toolbar7
    global figure7
    global ax7
    global ax6

    (w, h) = size
    inchsize = (w / 25.4, h / 25.4)
    figure7 = Figure(inchsize)

    ax6 = figure7.add_subplot(211, aspect='equal', xlabel='Re(Z)/\u03a9 cm²',
                             ylabel='-Im(Z) cm²')
    ax7 = figure7.add_subplot(212, xscale='log', xlabel='log(tau/s)',
                             ylabel='res')
    ax7.set_ylabel("h(tau)/Ohm cm²")


    # create canvas as matplotlib drawing area
    canvas_fit = FigureCanvasTkAgg(figure7, master=masterframe)
    toolbar7 = NavigationToolbar2Tk(canvas_fit, masterframe)
    canvas_fit.get_tk_widget().pack_forget()
    canvas_fit.get_tk_widget().pack()
    figure7.tight_layout()

def import_settings():
    global entry_header
    global clickede1
    global clickede2
    global clickede3
    global top01
    global collist
    global sheet
    global delimvalue



    top01 = Toplevel()
    top01.title('Import settings')
    top01.geometry("600x400+50+50")
    frame011 = Frame(top01)
    frame011.grid(row=0, column=0, padx=10, pady=10)

    presets_label = Label(frame011, text="Presets")
    presets_droplist = OptionMenu(frame011, presetvalue, *presetlist, command=update_importsetting)

    header_label = Label(frame011, text="Header lines")
    entry_header = Entry(frame011)
    entry_header.insert(0, str(headerv.get()))

    separator_label = Label(frame011, text="Delimiter")
    drop_separator = OptionMenu(frame011, delimvalue, *delimlist)

    drop1 = OptionMenu(frame011, clickede1, *collist)
    drop2 = OptionMenu(frame011, clickede2, *collist)
    drop3 = OptionMenu(frame011, clickede3, *collist)


    droplabel1 = Label(frame011, text="Frequency")
    droplabel2 = Label(frame011, text="Real part")
    droplabel3 = Label(frame011, text="Imaginary part")
    savebutton0 = tk.Button(frame011, text="Reload data", command=save_importsettings)
    closebutton = tk.Button(frame011, text="Close this window", command=top01.destroy)
    sheet = tksheet.Sheet(top01)

    #sheet.set_sheet_data(datav_list_all[0])

    presets_label.grid(row=0, column=0, padx=10, pady=10, sticky=NW)
    presets_droplist.grid(row=0, column=1, padx=10, pady=10, sticky=NW)
    header_label.grid(row=1, column=0, padx=10, pady=10, sticky=NW)
    entry_header.grid(row=1, column=1, padx=10, pady=10, sticky=NW)

    separator_label.grid(row=1, column=2, padx=10, pady=10, sticky=NW)
    drop_separator.grid(row=1, column=3, padx=10, pady=10, sticky=NW)

    droplabel1.grid(row=2, column=0, padx=10, pady=3)
    droplabel2.grid(row=2, column=1, padx=10, pady=3)
    droplabel3.grid(row=2, column=2, padx=10, pady=3)
    drop1.grid(row=3, column=0, padx=10, pady=3)
    drop2.grid(row=3, column=1, padx=10, pady=3)
    drop3.grid(row=3, column=2, padx=10, pady=3)

    savebutton0.grid(row=4, column=0, padx=10, pady=20, sticky=NW)
    closebutton.grid(row=4, column=1, padx=10, pady=20, sticky=NW)
    sheet.grid(row=5, column=0, padx=10, pady=20, sticky=NW)

def update_importsetting(event):
    if presetvalue.get() == "Zahner":
        clickede1.set(collist[1])
        clickede2.set(collist[2])
        clickede3.set(collist[3])
        entry_header.delete(0, END)
        entry_header.insert(0, str(12))
        headerv.set(12)
        delimvalue.set(delimlist[2])
    if presetvalue.get() == "Autolab":
        clickede1.set(collist[0])
        clickede2.set(collist[4])
        clickede3.set(collist[5])
        entry_header.delete(0, END)
        entry_header.insert(0, str(11))
        headerv.set(11)
        delimvalue.set(delimlist[0])
    if presetvalue.get() == "Biologic":
        clickede1.set(collist[2])
        clickede2.set(collist[0])
        clickede3.set(collist[1])
        entry_header.delete(0, END)
        entry_header.insert(0, str(0))
        headerv.set(0)
        delimvalue.set(delimlist[2])
    if presetvalue.get() == "Circuit Simulator":
        clickede1.set(collist[0])
        clickede2.set(collist[1])
        clickede3.set(collist[2])
        entry_header.delete(0, END)
        entry_header.insert(0, str(8))
        headerv.set(8)
        delimvalue.set(delimlist[2])

def save_importsettings():
    global clickede1
    global clickede2
    global clickede3
    global delimvalue
    global delimv
    global headerv
    global freqcol
    global realcol
    global imagcol
    global top01

    if entry_header.get() == 0:
        headerv.set("None")
    else:
        headerv.set(int(entry_header.get()))



    clickedev1 = clickede1.get()
    clickedev2 = clickede2.get()
    clickedev3 = clickede3.get()

    freqcol = collist.index(clickedev1)
    realcol = collist.index(clickedev2)
    imagcol = collist.index(clickedev3)

    #print((freqcol,realcol,imagcol))
    #top01.destroy()

    global text_list
    global datav_list_all
    global origin_list
    global df_list
    global instlist
    global labellist_all
    global delimiterv
    global sheet
    if not import_tree2.get_children():
        pass
    else:
        import_tree2.delete(*import_tree2.get_children())
    #text_files = filedialog.askopenfilenames()
    #text_list = list(text_files)

    count = 0
    labellist_all = []
    for record in text_list:
        path = pathlib.PurePath(record)
        import_tree2.insert(parent='', index='end', iid=count, text=path.name, values=record)
        labellist_all.append(path.name[:-4])
        count += 1

    origin_list = []
    x = import_tree2.get_children()
    child_id = import_tree2.get_children()[0]
    import_tree2.focus(child_id)
    import_tree2.selection_set(child_id)

    item_text = [text_list[i] for i in range(len(text_list))]
    if item_text[0][0][-4:] == '.mpr':
        files_all = [BioLogic.MPRfile(file) for file in item_text]
        df_all = [pd.DataFrame(files_all[i].data) for i in range(len(files_all))]
        datav_list_all = [data_list_i.values.tolist()[0:] for data_list_i in df_all]
    else:
        delimv = delimvalue.get()
        headv = headerv.get()
        data_list = [pd.read_csv(file, sep=delimv, header=headv, encoding='latin1') for file in item_text]
        datav_list_all = [data_list_i.values.tolist()[0:] for data_list_i in data_list]



    # import_settings()
    instlist = []
    try:
        create_instlists()
    except IndexError:
        #messagebox.showerror("Error", "Oops, something went wrong :(\nMaybe too few headlines?")
        top01.grab_set()
    import_tree2.bind('<ButtonRelease-1>', plotselection2)
    plotselection2(0)

    #sheet = tksheet.Sheet(top01)
    #sheet.grid(row=3, column=0)
    sheet.set_sheet_data(datav_list_all[0])

    # import_tree.bind('<ButtonRelease-1>', plotselection)
    #plotselection2(0)

def celldetails():
    global e_desc
    global e_area
    global e_thick
    global e_temp
    global sym_checkvar
    global top00
    global descv
    global areav
    global thickv
    global tempv
    global symcheckvarv
    global frame001


    top00 = Toplevel()
    top00.title('Cell details')
    top00.geometry("600x400+50+50")
    frame001 = Frame(top00)
    frame002 = Frame(top00)

    frame001.grid(row=0, column=0, padx=10, pady=10)
    frame002.grid(row=0, column=1, padx=10, pady=10)

    info_text = Label(frame002, text="Info:\nIn symmetrical cells, the normalized \nimpedance is calculated by dividing the \nmeasured impedance by 2 and multiplying \nit with the active area of one electrode.", justify=LEFT
          )
    info_text.grid(row=0, column=0, padx=10, pady=10)
    l_desc = Label(frame001, text="Description")
    e_desc = Entry(frame001)
    if descv:
        e_desc.insert(0, str(descv))


    sym_checkvar = IntVar()
    if symcheckvarv == 1:
        sym_checkvar.set(1)
    elif symcheckvarv == 0:
        sym_checkvar.set(0)
    l_symm = Label(frame001, text="Symmetrical cell?")
    symm_check = Checkbutton(frame001, variable=sym_checkvar)

    l_area = Label(frame001, text="Electrode area/cm²")
    e_area = Entry(frame001)
    e_area.insert(0, 1)
    if areav:
        e_area.delete(0, END)
        e_area.insert(0, str(areav))

    l_thick = Label(frame001, text="Separator thickness/µm")
    e_thick = Entry(frame001)
    if thickv:
        e_thick.insert(0, str(thickv))

    l_temp = Label(frame001, text="Temperature/°C")
    e_temp = Entry(frame001)
    if tempv:
        e_temp.insert(0, str(tempv))

    savebutton= tk.Button(frame001, text="Save and Close", command=savedetails)


    l_desc.grid(row=0, column=0, padx=3, pady=3)
    e_desc.grid(row=0, column=1, padx=3, pady=3)
    l_area.grid(row=1, column=0, padx=3, pady=3)
    e_area.grid(row=1, column=1, padx=3, pady=3)
    l_thick.grid(row=2, column=0, padx=3, pady=3)
    e_thick.grid(row=2, column=1, padx=3, pady=3)
    l_temp.grid(row=3, column=0, padx=3, pady=3)
    e_temp.grid(row=3, column=1, padx=3, pady=3)
    l_symm.grid(row=4, column=0, padx=3, pady=3)
    symm_check.grid(row=4, column=1, padx=3, pady=3)
    savebutton.grid(row=5, column=0, columnspan=2, padx=3, pady=3)

def savedetails():
    global descv
    global areav
    global thickv
    global tempv
    global symcheckvarv
    descv = e_desc.get()
    areav = e_area.get()
    thickv = e_thick.get()
    tempv = e_temp.get()
    symcheckvarv = sym_checkvar.get()
    top00.destroy()

def create_instlists():
    global instlist

    for i in range(len(datav_list_all)):
        instlist.append(DRTsolver(datav_list_all[i]))

def create_instlists2():
    global instlist2

    for i in range(len(datav_list_all2)):
        instlist2.append(Peakfit(datav_list_all2[i]))

def plotselection2(event):
    #delete_all_elements()
    ax1.cla()
    ax2.cla()

    x = import_tree2.selection()
    datav_list = [datav_list_all[int(item_i)] for item_i in x]
    for i in x:
        #getfminmax()
        #print(instlist[int(i)].freq)
        instlist[int(i)].plotdrt()

    plotresiduals()
    #getfminmax()

def plotselection3(event):
    x = import_drt_tree.selection()
    datav_list = [datav_list_all2[int(item_i)] for item_i in x]
    pw.clearplot()
    for i in x:
        #getfminmax()
        rtimesv = instlist2[int(i)].rtimes
        ampsn = instlist2[int(i)].ampsn
        realsv = instlist2[int(i)].realsn
        negimsv = instlist2[int(i)].negimn
        realsrec = instlist2[int(i)].reconrealsn
        negimsrec = instlist2[int(i)].reconnegimagsn
        pw.plotdrt(rtimesv,ampsn)
        pw.plotnyq(realsv,np.negative(negimsv))
        pw.plotrecon(realsrec, negimsrec)

        #print(negimsv)
        #print(realsv)


    cid = pw.canvas.mpl_connect('button_press_event', pw.onclick)
    pw.clearallpoints()
    points = pw.points
    pw.canvas.draw()

    #plotresiduals()

def savetxt():
    global file_name
    name_input = " "
    files = [('Text Document', '*.txt')]
    #file_name = filedialog.asksaveasfilename(filetypes=files, defaultextension=files, initialfile=name_input)
    file_name = filedialog.askdirectory()
    x = import_tree2.selection()
    for i in x:
        instlist[int(i)].exporttxt(labellist_all[int(i)])

def savetxt_area(exportlistv):
    global file_name2
    name_input = " "
    files = [('Text Document', '*.txt')]
    header = "Name\tR_serial/Ohm\tR_interface/Ohm\tR_diff/Ohm\tR_total/Ohm\n"
    file_name2 = filedialog.asksaveasfilename(filetypes=files, defaultextension=files, initialfile=name_input)

    #exportlistvv = np.column_stack(exportlistv)
    with open(file_name2, 'w') as f:
        f.write(header)
        for i in range(len(exportlistv)):
            f.write(exportlistv[i][0]+"\t"+exportlistv[i][1]+"\t"+exportlistv[i][2]+"\t"+exportlistv[i][3]+"\t"+exportlistv[i][4]+'\n')
    #np.savetxt('test.txt', ab, fmt="%10s %10.3f")
    #np.savetxt(file_name2[:-5] + "/" + "_areas" + ".txt", exportlistvv, header=header, delimiter='\t', fmt='%.5e')

def takesecond(elem):
    return elem[0]

def get_areas():
    x = import_drt_tree.selection()
    datav_list = [datav_list_all2[int(item_i)] for item_i in x]
    for i in x:
        rserial = instlist2[int(i)].rserial
        areas = pw.peakareas
        print(rserial,areas)

def quickfit():
    x = import_drt_tree.selection()
    datav_list = [datav_list_all2[int(item_i)] for item_i in x]
    pw.clearplot()
    logtau_range = slider_range.getValues()
    lowerblogv = logtau_range[0]
    middleblogv = logtau_range[1]
    for i in x:
        # getfminmax()
        rtimesv = instlist2[int(i)].rtimes
        logrtimes = instlist2[int(i)].logrtimes
        #print(logrtimes[-1])
        ampsn = instlist2[int(i)].ampsn
        amps = instlist2[int(i)].amps
        ampsn_pol = [1/(sum(amps))*ampsi for ampsi in amps]
        realsv = instlist2[int(i)].realsn
        negimsv = instlist2[int(i)].negimn
        realsrec = instlist2[int(i)].reconrealsn
        negimsrec = instlist2[int(i)].reconnegimagsn
        rserial = instlist2[int(i)].rserialn
        deltau = instlist2[int(i)].deltau
        pw.plotdrt(rtimesv, ampsn)
        pw.plotnyq(realsv, np.negative(negimsv))
        pw.plotrecon(realsrec, negimsrec)
        #print(rserial)
        #fitpeaks = instlist2[int(i)].fitpeaks
        lowerboundary = findclosest(logrtimes, lowerblogv)
        higherboundary = findclosest(logrtimes, middleblogv)
        pw.plotrectangle(lowerblogv, middleblogv)
        areaundercurve = AreaUnderCurve(logrtimes[lowerboundary:higherboundary], ampsn[lowerboundary:higherboundary])
        print(areaundercurve.calculate_area())
        #print(pw.pointelist)
        #fit1 = Peakfit(logrtimes, ampsn, pw.pointelist)
        #pw.plotfit()
    pw.canvas.draw()
        #areas = pw.printpeakareas()
        #print(areas)
        #print(deltau)
        #totalresistance = sum(amps) + rserial[0]
        #areas2 = [areas_i/deltau[0]*totalresistance for areas_i in areas]
        #print(areas2)

        #print(rserial)
        #pw.peakareas()
        #pw.get_area(rserial)
        #pw.plotnyq(realsv, np.negative(negimsv))
        #pw.plotrecon(realsrec, negimsrec)
    #get_areas()

def areafit():
    x = import_drt_tree.selection()
    datav_list = [datav_list_all2[int(item_i)] for item_i in x]

    global labellist_all2
    # print(labellist_all)

    pw.clearplot()

    label_list = []
    rserial_list = []
    rinterface_list = []
    rdiff_list = []
    rtotal_list= []
    for i in x:
        rtimesv = instlist2[int(i)].rtimes
        logrtimes = instlist2[int(i)].logrtimes
        lowerboundary_index = findclosest(logrtimes,lowerblog)
        middleboundary_index = findclosest(logrtimes,middleblog)

        ampsn = instlist2[int(i)].ampsn
        amps = instlist2[int(i)].amps
        rserialn = instlist2[int(i)].rserialn[0]
        deltau = instlist2[int(i)].deltau

        right_area = sum(ampsn[middleboundary_index:])
        middle_area = sum(ampsn[lowerboundary_index:middleboundary_index])
        left_area = sum(ampsn[:lowerboundary_index])

        label_list.append(str(labellist_all2[int(i)]))
        rserial_list.append(left_area+rserialn)
        rinterface_list.append(middle_area)
        rdiff_list.append(right_area)
        rtotal_list.append(right_area+left_area+rserialn+middle_area)

        #print(right_area,left_area+rserial)
    newexport_table = [[label_list[i],str(rserial_list[i]), str(rinterface_list[i]),str(rdiff_list[i]),str(rtotal_list[i])] for i in range(len(label_list))]
    #print(newexport_table)
    #exportlist = [label_list,rserial_list,rinterface_list,rtotal_list]
    print(newexport_table)
    #savetxt_area(newexport_table)

def popzero(liste):
    #listeneu = [i for i in liste if i != 0]
    listeneu = liste
    return listeneu

def findclosest(listv,x):
    index = min(range(len(listv)), key=lambda i: abs(listv[i] - x))
    return index

def areaconfirm():
    global arean
    arean = float(area_entry.get())
    plotselection2(0)

class Peakfit:
    def __init__(self, data):
        self.data = data
        self.model = None
        # initialize original freq, reals, negims
        self.freq = popzero([col[0] for col in self.data])
        # reverse freqlist, if it starts from low to high values

        if np.mean([col[6] for col in self.data]) < 0:
            self.negimn = popzero([col[6] for col in self.data])
        else:
            self.negimn = popzero([-col[6] for col in self.data])
        self.realsn = popzero([col[5] for col in self.data])
        # print(self.reals)
        self.reconrealsn = popzero([col[7] for col in self.data])
        self.reconnegimagsn = popzero([col[8] for col in self.data])
        self.rtimes = [col[11] for col in self.data]
        self.logrtimes = [col[12] for col in self.data]
        self.amps = [col[13] for col in self.data]
        self.ampsn = [col[14] for col in self.data]
        self.rserialn = [[col[16] for col in self.data][0]]
        # print(self.rserial)
        self.cserial = popzero([col[17] for col in self.data])
        self.lserial = popzero([col[18] for col in self.data])
        self.deltau = popzero([col[21] for col in self.data])
        self.elarea = popzero([col[22] for col in self.data])

    def create_multipeak_model(self, num_peaks):
        model = GaussianModel(prefix='p0_')
        for i in range(1, num_peaks):
            model += GaussianModel(prefix=f'p{i}_')
        return model

    def fit_spectrum(self):
        if self.data is not None:
            x_data_fit = np.linspace(min(self.logrtimes),max(self.logrtimes),1000)
            self.model = self.create_multipeak_model(3)
            params = self.model.make_params()


    def fit_fct(self):
        params = self.model.make_params()
        for i, (x0, y0) in enumerate(self.points):
            params[f'p{i}_center'].set(value=x0, min=self.xdata[0], max=self.xdata[-1])
            params[f'p{i}_sigma'].set(value=0.1, min=0.01, max=0.5)
            params[f'p{i}_amplitude'].set(value=y0, min=0)
        result = self.model.fit(self.ydata, params, x=self.xdata)
        return result

    def fitvalues(self):
        yfit = self.fit_fct.best_values
        bestresult = (self.xdata, yfit)
        return bestresult

    def peakareas(self):
        peak_areas = []
        for i in range(len(self.points)):
            area = self.fit_fct.best_values[f'p{i}_amplitude'] * np.sqrt(2 * np.pi) * self.fit_fct.best_values[f'p{i}_sigma']
            peak_areas.append(area)

    def individual_peaks(self):
        peakfctlist = []
        for i in range(len(self.points)):
            result = self.fit_fct
            yfit = result.eval_components()[f'p{i}_']
            peakfctlist.append([self.xdata, yfit])
        return peakfctlist

class AreaUnderCurve:
    def __init__(self, x_data, y_data):
        self.x_data = x_data
        self.y_data = y_data
        self.num_points = len(x_data)

    def calculate_area(self):
        area = 0
        for i in range(self.num_points - 1):
            area += (self.y_data[i] + self.y_data[i + 1]) * (self.x_data[i + 1] - self.x_data[i]) / 2
        return area

class Plotwindow:
    def __init__(self, masterframe, size):
        self.size = size
        self.masterframe = masterframe

        (w, h) = size
        inchsize = (w / 25.4, h / 25.4)
        self.figure = Figure(inchsize)
        self.axes1 = self.figure.add_subplot(211, aspect='equal', xlabel='Re(Z)/\u03a9 cm²',
                                  ylabel='-Im(Z)/\u03a9 cm²')
        self.axes2 = self.figure.add_subplot(212, xscale='log', xlabel=r'log($\tau$/s)',
                                  ylabel='resist')

        self.axes2.set_ylabel(r"$h(\tau)/\mathrm{\Omega cm^2}$")

        self.axes3 = self.axes2.secondary_xaxis('top')
        self.axes3.set_xlabel('Upper X-axis')

        # create canvas as matplotlib drawing area
        self.canvas = FigureCanvasTkAgg(self.figure, master=masterframe)
        self.canvas.get_tk_widget().pack()

        self.toolbar = NavigationToolbar2Tk(self.canvas, masterframe)

        self.points = []
        self.pointelist = []
        self.old_coords = None
        # self.line1 = self.createline(0,0,0,0)
        self.cid0 = self.canvas.mpl_connect('button_press_event', self.onclick)
        self.peakareas = []
        self.vspan = self.axes2.axvspan(0, 0, alpha=0.5, color='grey')
        # self.blinepoints = []

    def plotdrt(self, x, y):

        self.axes2.plot(x, y)
        self.axes2.set_ylabel(r"$h(\tau)/\mathrm{\Omega cm^2}$")
        self.axes2.set_xlabel(r"$\tau/\mathrm{s}$")
        self.axes2.set_xscale('log')
        self.axes3 = self.axes2.secondary_xaxis('top')
        self.axes3.set_xlabel('Upper X-axis')
        self.figure.tight_layout()
        self.canvas.draw()

    def plotnyq(self,x,y):
        self.figure.tight_layout()
        color = next(self.axes1._get_lines.prop_cycler)['color']
        self.axes1.set_ylabel('Re(Z)/\u03a9 cm²')
        self.axes1.set_xlabel('-Im(Z)/\u03a9 cm²')
        self.axes1.plot(x,y, "o", markersize=4, markerfacecolor='none', color=color)
        self.canvas.draw()
    def plotrecon(self,x,y):
        color = next(self.axes1._get_lines.prop_cycler)['color']
        self.axes1.plot(x, y, "-", markersize=4, markerfacecolor='none',
                 color=color)
        self.canvas.draw()
    def plotrectangle(self,lowerboundary,upperboundary):
        self.vspan.remove()
        self.vspan = self.axes2.axvspan(10**lowerboundary, 10**upperboundary, alpha=0.5, color='grey')
        self.canvas.draw()
    def clearplot(self):
        self.axes1.cla()
        self.axes2.cla()
        self.canvas.draw()
    def onclick(self, event):
        if event.button == 1:
            x, y = event.xdata, event.ydata
            if x is not None and y is not None:
                self.pointelist.append((x,y))
                self.axes2.plot(x, y, 'ro', markersize=8)
                self.canvas.draw()

    def clearpoints(self):
        self.pointelist[-1].remove()
        del self.pointelist[-1]
        del self.points[-1]
        self.canvas.draw()

    def clearallpoints(self):
        self.pointelist = []
        self.canvas.draw()

    def printpoints(self):
        print(self.pointelist)


    def firstclick(self, event):
        #if str(event.type) == 'button_press_event':
        #x,y = event.xdata, event.ydata
        #self.line1.remove()
        self.old_coords = [event.xdata, event.ydata]
        #self.line1.remove()
        self.canvas.draw()
        self.cid0 = self.canvas.mpl_connect('motion_notify_event', self.drawfct)

    def secondclick(self, event):
        x, y = event.xdata, event.ydata
        x1, y1 = self.old_coords
        #print(x,y)
        #print(x1,y1)
        self.canvas.mpl_disconnect(self.cid0)

    def drawfct(self, event):
        self.blinepoints = []
        self.line1.remove()
        x, y = event.xdata, event.ydata
        if self.old_coords:
            x1, y1 = self.old_coords
            self.line1 = self.createline(x, y, x1, y1)
            self.blinepoints.append(x)
            self.blinepoints.append(y)
            self.blinepoints.append(x1)
            self.blinepoints.append(y1)

class DRTsolver:
    def __init__(self, data):
        self.data = data
        # initialize original freq, reals, negims
        self.freqi = [col[freqcol] for col in self.data]
        # reverse freqlist, if it starts from low to high values
        if np.mean([col[imagcol] for col in self.data]) < 0:
            self.negimi = [col[imagcol] for col in self.data]
        else:
            self.negimi = [-col[imagcol] for col in self.data]
        self.realsi = [col[realcol] for col in self.data]

        if self.freqi[0] < self.freqi[1]:
            self.freqi.reverse()
            self.realsi.reverse()
            self.negimi.reverse()
        # initialize rfact and lambada
        self.rfact = 2
        self.lambada = 0.2
        self.trufals = True
        self.hboundary = 0
        self.lboundary = len(self.freqi)
        self.boundaryfactor = 10

        # initialize adjusted freq, negim and reals
        self.freq = self.freqi[self.hboundary:self.lboundary]
        self.negim = self.negimi[self.hboundary:self.lboundary]
        self.reals = self.realsi[self.hboundary:self.lboundary]

    def findclosest(self, lvalue, hvalue):
        closestlvalue = min(self.freqi, key=lambda x: abs(x - lvalue))
        closesthvalue = min(self.freqi, key=lambda x: abs(x - hvalue))
        self.hboundary = self.freqi.index(closesthvalue)
        self.lboundary = self.freqi.index(closestlvalue)

        self.freq = self.freqi[self.hboundary:self.lboundary]
        self.negim = self.negimi[self.hboundary:self.lboundary]
        self.reals = self.realsi[self.hboundary:self.lboundary]
        #print(max(self.freq), min(self.freq))
    def findclosesth(self, value):
        closestvalue = min(self.freq, key=lambda x: abs(x - value))
        self.hboundary = self.freq.index(closestvalue)
        self.freq = self.freqi[self.hboundary:self.lboundary]
        self.negim = self.negimi[self.hboundary:self.lboundary]
        self.reals = self.realsi[self.hboundary:self.lboundary]

    def findclosestl(self, value):
        closestvalue = min(self.freq, key=lambda x: abs(x - value))
        self.lboundary = self.freq.index(closestvalue)
        self.freq = self.freqi[self.hboundary:self.lboundary]
        self.negim = self.negimi[self.hboundary:self.lboundary]
        self.reals = self.realsi[self.hboundary:self.lboundary]

    def lowerhboundary(self):
        if self.hboundary < len(self.data):
            self.hboundary = self.hboundary + 1
            #print("h=" + str(self.hboundary))
            self.freq = self.freqi[self.hboundary:self.lboundary]
            self.negim = self.negimi[self.hboundary:self.lboundary]
            self.reals = self.realsi[self.hboundary:self.lboundary]

    def higherhboundary(self):
        if self.hboundary > 0:
            self.hboundary = self.hboundary - 1
            #print("h="+ str(self.hboundary))
            self.freq = self.freqi[self.hboundary:self.lboundary]
            self.negim = self.negimi[self.hboundary:self.lboundary]
            self.reals = self.realsi[self.hboundary:self.lboundary]

    def higherlboundary(self):
        if self.lboundary > 0:
            self.lboundary = self.lboundary - 1
            #print("l=" + str(self.lboundary))
            self.freq = self.freqi[self.hboundary:self.lboundary]
            self.negim = self.negimi[self.hboundary:self.lboundary]
            self.reals = self.realsi[self.hboundary:self.lboundary]

    def lowerlboundary(self):
        if self.lboundary < len(self.data):
            self.lboundary = self.lboundary + 1
            #print("l=" + str(self.lboundary))
            self.freq = self.freqi[self.hboundary:self.lboundary]
            self.negim = self.negimi[self.hboundary:self.lboundary]
            self.reals = self.realsi[self.hboundary:self.lboundary]

    def fminmax(self):
        fmin = min(self.freq)
        fmax = max(self.freq)
        return fmin,fmax
    def get_rfact(self):
        if drop.get() == options[0]:
            self.rfact = 1
        elif drop.get() == options[1]:
            self.rfact = 2
        elif drop.get() == options[2]:
            self.rfact = 3
        elif drop.get() == options[3]:
            self.rfact = 4
        elif drop.get() == options[4]:
            self.rfact = 6
        else:
            self.rfact = 1
    def get_boundaryfactor(self):
        if drop_b.get() == boundaryoptions[0]:
            self.boundaryfactor = 1
        if drop_b.get() == boundaryoptions[1]:
            self.boundaryfactor = 10
        if drop_b.get() == boundaryoptions[2]:
            self.boundaryfactor = 100
        if drop_b.get() == boundaryoptions[3]:
            self.boundaryfactor = 1000
    def get_lambda(self):
        self.lambada = float(e1.get())

    def rtimesfct(self):
        self.get_boundaryfactor()
        lowerboundary = (1/self.boundaryfactor)*1 / (2 * math.pi * self.freq[0])
        upperboundary = self.boundaryfactor*1 / (2 * math.pi * self.freq[-1])
        self.get_rfact()
        rtimes = np.geomspace(lowerboundary, upperboundary, num=self.rfact * len(self.freq)).tolist()
        return rtimes

    def buildmat(self):
        if self.trufals:
            self.get_lambda()
        rtimes = self.rtimesfct()
        realpartlist = [[(1 / (1 + 2j * math.pi * i * j)).real for j in rtimes] for i in self.freq]
        imagpartlist = [[(1 / (1 + 2j * math.pi * i * j)).imag for j in rtimes] for i in self.freq]
        matA = np.vstack([realpartlist, imagpartlist]).tolist()
        vartp = [var_r.get(), var_c.get(), var_l.get()]
        regmat1 = [[x * self.lambada for x in y] for y in np.identity(len(list(rtimes))).tolist()]
        if vartp == [1, 1, 1]:
            r_array = [[1, 0, 0] for i in self.freq]
            c_array = [[0, -1 / (2 * math.pi * i), 2 * math.pi * i] for i in self.freq]
            rc_array = np.vstack([r_array, c_array]).tolist()
            matC = np.hstack([matA, rc_array]).tolist()
            regmat2 = [[0, 0, 0] for i in rtimes]
        elif vartp == [1, 0, 0]:
            r_array = [[1] for i in self.freq]
            c_array = [[0] for i in self.freq]
            rc_array = np.vstack([r_array, c_array]).tolist()
            matC = np.hstack([matA, rc_array]).tolist()
            regmat2 = [[0] for i in rtimes]
        elif vartp == [1, 1, 0]:
            r_array = [[1, 0] for i in self.freq]
            c_array = [[0, -1 / (2 * math.pi * i)] for i in self.freq]
            rc_array = np.vstack([r_array, c_array]).tolist()
            matC = np.hstack([matA, rc_array]).tolist()
            regmat2 = [[0, 0] for i in rtimes]
        elif vartp == [0, 1, 0]:
            r_array = [[0] for i in self.freq]
            c_array = [[-1 / (2 * math.pi * i)] for i in self.freq]
            rc_array = np.vstack([r_array, c_array]).tolist()
            matC = np.hstack([matA, rc_array]).tolist()
            regmat2 = [[0] for i in rtimes]
        elif vartp == [0, 0, 1]:
            r_array = [[0] for i in self.freq]
            c_array = [[2 * math.pi * i] for i in self.freq]
            rc_array = np.vstack([r_array, c_array]).tolist()
            matC = np.hstack([matA, rc_array]).tolist()
            regmat2 = [[0] for i in rtimes]
        elif vartp == [1, 0, 1]:
            r_array = [[1, 0] for i in self.freq]
            c_array = [[0, 2 * math.pi * i] for i in self.freq]
            rc_array = np.vstack([r_array, c_array]).tolist()
            matC = np.hstack([matA, rc_array]).tolist()
            regmat2 = [[0, 0] for i in rtimes]
        elif vartp == [0, 1, 1]:
            r_array = [[0, 0] for i in self.freq]
            c_array = [[-1 / (2 * math.pi * i), 2 * math.pi * i] for i in self.freq]
            rc_array = np.vstack([r_array, c_array]).tolist()
            matC = np.hstack([matA, rc_array]).tolist()
            regmat2 = [[0, 0] for i in rtimes]
        else:
            r_array = [[0] for i in self.freq]
            c_array = [[0] for i in self.freq]
            rc_array = np.vstack([r_array, c_array]).tolist()
            matC = np.hstack([matA, rc_array]).tolist()
            regmat2 = [[0] for i in rtimes]
        regmat3 = np.hstack([regmat1, regmat2]).tolist()
        matCreg = np.vstack([matC, regmat3])
        return matCreg

    def buildbvec(self):
        rtimes = self.rtimesfct()
        bvec = np.hstack([self.reals, self.negim]).tolist()
        regvec0 = [0 for i in rtimes]
        bvecreg = np.hstack([bvec, regvec0])
        return bvecreg

    def nnlsfunc(self):
        mat = self.buildmat()
        vec = self.buildbvec()
        result = nnls(mat, vec)[0].tolist()
        return result

    def multilambda(self):
        #print(lamalist)
        listofpoints = []
        self.trufals = False
        for itemi in lamalist:
            self.lambada = itemi
            listofpoints.append(self.residuals())
        self.trufals = True
        return listofpoints

    def residuals(self):
        mat = self.buildmat()
        vec = self.buildbvec()
        resnorm = nnls(mat, vec)[1]
        solnorm = np.linalg.norm(self.amplitudes())
        norms = [resnorm, solnorm]
        return norms

    def plotres(self):
        #ax4.set_xscale('log')
        #ax5.set_xscale('log')
        ax4.set_xscale('linear')
        ax4.cla()
        ax5.cla()
        ax4.set_ylabel("real residual /%")
        ax5.set_ylabel("imaginary residual /%", color='C1')
        ax5.yaxis.set_label_position('right')
        residu = self.residualk()
        middle_item = len(residu) // 2
        realresidu = self.residualk()[:middle_item]
        negimresidu = self.residualk()[middle_item:]

        ax4.plot(self.freq, realresidu, "-o", markersize=4, markerfacecolor='none')
        ax5.plot(self.freq, negimresidu, "-o", markersize=4, color="C1", markerfacecolor='none')
        ax4.set_xscale('log')
        ax5.set_xscale('log')

        ax5.tick_params(axis='y', colors='C1', which='both')
        ax4.format_coord = lambda x, y: ''
        ax5.format_coord = lambda x, y: ''


        ax5.spines['right'].set_color('C1')
        ax4.set_xlabel("Frequency / Hz")
        figurek.tight_layout()
        canvas4.draw()

    def rsquared(self):
        mat = self.buildmat()
        vec = self.buildbvec()
        resnorm = nnls(mat, vec)[1]
        meanvalue = np.mean(self.amplitudes())
        sumofsquares = np.sum([(i-meanvalue)**2 for i in self.amplitudes()])
        rsquared = 1- (resnorm/sumofsquares)
        return rsquared

    def amplitudes(self):
        result = self.nnlsfunc()
        vartp = [var_r.get(), var_c.get(), var_l.get()]
        svartp = sum(vartp)
        if svartp == 0:
            n_result = result[:-1]
        else:
            n_result = result[:-svartp]
        return n_result

    def plotdrt(self):
        #ax1.cla()
        #ax2.cla()
        x = import_tree2.selection()
        labellist = []
        global labellist_all

        for xi in x:
            labellist.append(str(labellist_all[int(xi)]))
        #labellist.append("test")
        # set color of points and fits equal
        color = next(ax1._get_lines.prop_cycler)['color']
        realss=[arean*itemx for itemx in self.reals]
        realssrec = [arean * itemx for itemx in self.reconstruct()[0]]
        negimagss=[arean*itemx for itemx in np.negative(self.negim)]
        negimagssrec=[arean*itemx for itemx in np.negative(self.reconstruct()[1])]
        ax1.plot(realss, negimagss, "o", markersize=4, markerfacecolor='none', color= color)
        ax1.plot(realssrec, negimagssrec, "-", markersize=4, markerfacecolor='none', color = color)

        #inst1 = DRTsolver(freq, reals, negim, rfac, lama)
        ax1.set_xlabel('Re(Z)/\u03a9 cm²')
        ax1.set_ylabel('-Im(Z)/\u03a9 cm²')
        n_result = self.amplitudes()
        rtimes = self.rtimesfct()

        #print(sum(n_result))
        self.show_rcl()
        #print(self.rsquared())
        #print(self.get_rcl())
        delta_r = 1/(abs(np.log10(rtimes[10]))-abs(np.log10(rtimes[11])))
        #delta_r2 = 1/(abs(np.log(rtimes[100]))-abs(np.log(rtimes[101])))

        #result = self.nnlsfunc()
        y_data_n = [item / sum(n_result) for item in n_result]
        y_data = [arean*delta_r*item for item in n_result]
        if var_n.get() == 1:
            ax2.plot(rtimes, y_data_n)
        else:
            ax2.plot(rtimes, y_data)
        ax2.set_xscale('log')
        ax2.set_ylabel('h/\u03a9 cm²')
        ax2.set_xlabel(r'$ \tau$/s')
        ax2t.set_xlabel('f/Hz')
        ax2t.set_xscale('log')
        ax2t.xaxis.set_major_formatter(
            FuncFormatter(lambda y, _: '$10^{{{}}}$'.format(int(np.log10(1 / 2 * np.pi * y)))))
        ax2.set_ylabel("h(tau)/Ohm cm2")
        if var_legend.get() == 0:
            if len(labellist) < 10:
                ax2.legend(labellist)
        figure.tight_layout()
        self.reconstruct()
        canvas2.draw()

    def reconstruct(self):
        rtimeslength = len(np.array(self.rtimesfct()))
        freqlength = len(self.reals)
        bmatresult = self.buildmat()
        zdrt = np.array(bmatresult).dot(np.array(self.nnlsfunc()))[:2*freqlength]
        zmeas = self.buildbvec()[:2*freqlength]
        realsdrt = zdrt[:freqlength]
        imagsdrt = zdrt[freqlength:]
        return [realsdrt,imagsdrt]

    def residualk(self):
        rtimeslength = len(np.array(self.rtimesfct()))
        freqlength = len(self.reals)
        bmatresult = self.buildmat()
        zdrt = np.array(bmatresult).dot(np.array(self.nnlsfunc()))[:2 * freqlength]
        zmeas = self.buildbvec()[:2 * freqlength]
        realsdrt = zdrt[:freqlength]
        imagsdrt = zdrt[freqlength:]

        middle_item = len(list(zmeas)) // 2

        zabs0 = [np.sqrt(xx ** 2 + yy ** 2) for xx, yy in zip(zmeas[middle_item:], zmeas[:middle_item])]
        zabs1 = zabs0 + zabs0
        residu = np.array([100 * (np.array(zmeas)[i] - zdrt[i]) / np.array(zabs1)[i] for i in range(len(zdrt))])
        return residu

    def show_rcl(self):
        result = self.nnlsfunc()
        vartp = [var_r.get(), var_c.get(), var_l.get()]
        svartp = sum(vartp)
        if vartp == [1, 1, 1]:

            if result[-svartp+1] > 0.000001:
                r_label.config(text="R=" + str(float('%.3g' % result[-svartp])) + " Ohm")
                c_label.config(text="C=" + str(float('%.3g' % (1/result[-svartp+1])))+ " F")
                l_label.config(text="L=" + str(float('%.3g' % result[-svartp + 2])) + " H")
                rcl_list = [float('%.3g' % result[-svartp]),
                            float('%.3g' % (1 / result[-svartp + 1])),
                            float('%.3g' % result[-svartp + 2])]
            else:
                c_label.config(text="C=NaN")
                l_label.config(text="L=" + str(float('%.3g' % result[-svartp+2]))+ " H")
                rcl_list = [float('%.3g' % result[-svartp]),
                            float("NaN"),
                            float('%.3g' % result[-svartp+2])]
        elif vartp == [1, 0, 0]:
            r_label.config(text="R=" + str(float('%.3g' %result[-svartp])) + " Ohm")
            c_label.config(text="C=" + str(0) + " F")
            l_label.config(text="L=" + str(0) + " H")
            rcl_list = [float('%.3g' %result[-svartp]),
                        0,
                        0]
        elif vartp == [1, 1, 0]:
            r_label.config(text="R=" + str(float('%.3g' %result[-svartp])) + " Ohm")
            c_label.config(text="C=" + str(float('%.3g' % (1 / result[-svartp + 1]))) + " F")
            l_label.config(text="L=" + str(0) + " H")
            rcl_list = [float('%.3g' %result[-svartp]),
                        float('%.3g' % (1 / result[-svartp + 1])),
                        0]
        elif vartp == [0, 1, 0]:
            r_label.config(text="R=" + str(0) + " Ohm")
            c_label.config(text="C=" + str(float('%.3g' % (1 / result[-svartp]))) + " F")
            l_label.config(text="L=" + str(0) + " H")
            rcl_list = [0,
                        float('%.3g' % (1 / result[-svartp])),
                        0
                        ]
        elif vartp == [0, 0, 1]:
            r_label.config(text="R=" + str(0) + " Ohm")
            c_label.config(text="C=" + str(0) + " F")
            l_label.config(text="L=" + str(float('%.3g' %result[-svartp])) + " H")
            rcl_list = [0,
                        0,
                        float('%.3g' % result[-svartp])
                        ]
        elif vartp == [1, 0, 1]:
            r_label.config(text="R=" + str(float('%.3g' %result[-svartp])) + " Ohm")
            c_label.config(text="C=" + str(0) + " F")
            l_label.config(text="L=" + str(float('%.3g' % result[-svartp + 1])) + " H")
            rcl_list = [float('%.3g' %result[-svartp]),
                        0,
                        float('%.3g' % result[-svartp + 1])
                        ]
        elif vartp == [0, 1, 1]:
            r_label.config(text="R=" + str(0) + " Ohm")
            c_label.config(text="C=" + str(float('%.3g' % (1 / result[-svartp]))) + " F")
            l_label.config(text="L=" + str(float('%.3g' % result[-svartp + 1])) + " H")
            rcl_list = [0,
                        float('%.3g' % (1 / result[-svartp])),
                        float('%.3g' % result[-svartp + 1])]
        else:
            r_label.config(text="R=" + str(0) + " Ohm")
            c_label.config(text="C=" + str(0) + " F")
            l_label.config(text="L=" + str(0) + " H")
            rcl_list=[0,0,0]
        return rcl_list

    def exporttxt(self,labelname):
        arean = float(area_entry.get())
        residu = self.residualk()
        middle_item = len(residu) // 2
        realresidu = self.residualk()[:middle_item]
        negimresidu = self.residualk()[middle_item:]
        n_resulte = self.amplitudes()
        rtimese = self.rtimesfct()
        delta_logr = abs(np.log10(rtimese[10])) - abs(np.log10(rtimese[11]))
        delta_inv = 1 / delta_logr
        logrtimes = [np.log10(x) for x in rtimese]
        #y_data = [item for item in n_resulte]
        #y_data_a = [arean*item for item in y_data]
        y_data_n = [delta_inv*item for item in n_resulte]
        y_data_na = [arean*item for item in y_data_n]
        rcl_list = self.show_rcl()

        realsn = [arean*item for item in self.reals]
        negimn = [arean*item for item in self.negim]
        recon_real_n = [arean*item for item in self.reconstruct()[0].tolist()]
        recon_negim_n = [arean*item for item in self.reconstruct()[1].tolist()]
        expdatlist = [self.freq,
                   self.reals,
                   np.negative(self.negim).tolist(),
                   self.reconstruct()[0].tolist(),
                   np.negative(self.reconstruct()[1]).tolist(),
                   realsn,
                   np.negative(negimn).tolist(),
                   recon_real_n,
                   np.negative(recon_negim_n).tolist(),
                   realresidu.tolist(),
                   negimresidu.tolist(),
                   rtimese,
                   logrtimes,
                   y_data_n,
                   y_data_na,
                   [rcl_list[0]],
                   [arean*rcl_list[0]],
                   [rcl_list[1]],
                   [rcl_list[2]],
                   [self.lambada],
                   [self.rfact * len(self.freq)],
                   [delta_logr],
                   [arean]
                   ]
        #print(expdatlist)
        def make_same_length(listoflists):
            n = 0
            lenlist=[]
            for y in listoflists:
                lenlist.append(len(y))
            newlist = []
            for x in listoflists:
               if len(x) < max(lenlist):
                    x = x + [0]*(max(lenlist)-len(x))
                    newlist.append(x)
               else:
                   newlist.append(x)
                    #x = [item for sublist in x for item in sublist]
            return newlist
        expdatlistn = make_same_length(expdatlist)
        #print(len(expdatlistn))
        #print(expdatlistn)
        header = ["f/Hz",
                  "Re(Z)/Ohm",
                  "-Im(Z)/Ohm",
                  "Re(Z)fit/Ohm",
                  "-Im(Z)fit/Ohm",
                  "Re_norm/Ohm cm2",
                  "-Im_norm/Ohm cm2",
                  "Re_fit_norm/Ohm cm2",
                  "-Im_fit_norm/Ohm cm2",
                  "Real residual/%",
                  "Imaginary residual/%",
                  "tau/s",
                  "log(tau/s)",
                  "h(tau)/Ohm",
                  "h(tau)/Ohm cm2",
                  "Rserial/Ohm",
                  "Rserial/Ohm cm2",
                  "Cserial/F",
                  "Lserial/H",
                  "lambda",
                  "N_tau",
                  "DeltaLogTau",
                  "Area/cm2"]

        dfexp = pd.DataFrame(expdatlist)
        dfexp = dfexp.transpose()
        dfexp.columns = header
        #print(dfexp.to_string())
        #dfexp.to_csv(file_name[:-5]+"/"+labelname+"_drt"+".txt", sep='\t', index=False)
        file_path = file_name + "/" + labelname + "_drt" + ".txt"
        if os.path.exists(file_path):
            overwrite = messagebox.askyesno("File exists", f"The file name '{labelname}_drt.txt' already exists. Do you want to overwrite it?")
            if not overwrite:
                return

        dfexp.to_csv(file_path, sep='\t', index=False)
        #print(dfexp)

def plotresiduals():
    x = import_tree2.selection()
    instlist[int(x[0])].plotres()

def plotlamaopti():
    #entryset = [0.001, 1, 30]
    #lamalist = np.geomspace(entryset[0], entryset[1], entryset[2])
    #norms = np.array(lamaopti(entryset[0], entryset[1], entryset[2]))
    x = import_tree2.selection()
    datav_list = [datav_list_all[int(item_i)] for item_i in x]
    pointslist = DRTsolver(datav_list[0]).multilambda()

    #print(pointslist)
    #print(norms)
    #rsquared = np.array(lamaoptir(entryset[0], entryset[1], entryset[2]))
    #print(norms)
    '''
    resnorm = norms[:,0]
    solnorm = norms[:,1]
    fig, ax = plt.subplots()
    dat, = plt.plot(solnorm, resnorm, 'o-', picker=8)
    annot = ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Residual norm")
    plt.ylabel("Solution norm")
    names = [str(round(lamalist[n], 3)) for n in range(len(lamalist))]

    def update_annot(ind):
        x, y = dat.get_data()
        annot.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
        text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = dat.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    def onpick1(event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            print('onpick1 line:', np.column_stack([xdata[ind], ydata[ind]]))

    # fig.canvas.mpl_connect('pick_event', onpick1)

    fig.canvas.mpl_connect("motion_notify_event", hover)

    plt.show()
    '''
    # except: messagebox.showerror("Error", "No data imported")

def lcurveplot():
    # try:
    toplambda = Toplevel()
    toplambda.title("Lambda range")
    x = root.winfo_x()
    y = root.winfo_y()
    toplambda.geometry("+%d+%d" % (x + 500, y + 500))
    toplambda.focus_force()
    #toplambda.geometry("300x200")
    label1 = Label(toplambda, text="Start value")
    label2 = Label(toplambda, text="End value")
    label3 = Label(toplambda, text="Number of points")
    label1.grid(row=0, column=0)
    label2.grid(row=0, column=1)
    label3.grid(row=0, column=2)
    lambdaentry1 = Entry(toplambda)
    lambdaentry2 = Entry(toplambda)
    lambdaentry3 = Entry(toplambda)
    lambdaentry1.insert(0, "0.001")
    lambdaentry2.insert(0, "10")
    lambdaentry3.insert(0, "30")
    lambdaentry1.grid(row=1, column=0)
    lambdaentry2.grid(row=1, column=1)
    lambdaentry3.grid(row=1, column=2)

    def okcommand():
        global lamalist
        lamalist = np.geomspace(float(lambdaentry1.get()), float(lambdaentry2.get()), int(lambdaentry3.get()))
        toplambda.destroy()

        x = import_tree2.selection()
        datav_list = [datav_list_all[int(item_i)] for item_i in x]
        pointslist = DRTsolver(datav_list[0]).multilambda()
        # fig = plt.figure(3, figsize=(6, 4), dpi=80)
        fig, ax = plt.subplots()
        resnorm = [i[0] for i in pointslist]
        solnorm = [i[1] for i in pointslist]
        dat, = plt.plot(resnorm, solnorm, 'o-', picker=8)
        annot = ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Residual norm")
        plt.ylabel("Solution norm")
        names = [str(round(lamalist[n], 3)) for n in range(len(lamalist))]

        def update_annot(ind):
            x, y = dat.get_data()
            annot.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
            text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
            annot.set_text(text)
            annot.get_bbox_patch().set_alpha(0.4)

        def hover(event):
            vis = annot.get_visible()
            if event.inaxes == ax:
                cont, ind = dat.contains(event)
                if cont:
                    update_annot(ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

        def onpick1(event):
            if isinstance(event.artist, Line2D):
                thisline = event.artist
                xdata = thisline.get_xdata()
                ydata = thisline.get_ydata()
                ind = event.ind
                print('onpick1 line:', np.column_stack([xdata[ind], ydata[ind]]))

        # fig.canvas.mpl_connect('pick_event', onpick1)

        fig.canvas.mpl_connect("motion_notify_event", hover)

        plt.show()

        # except: messagebox.showerror("Error", "No data imported")

    okbutton = Button(toplambda, text="OK", command=okcommand)
    okbutton.grid(row=2, column=1)

    cancelbutton = Button(toplambda, text="Cancel", command=toplambda.destroy)
    cancelbutton.grid(row=2, column=2)

def residualplot():
    # try:
    # try:
    toplambda = Toplevel()
    toplambda.title("Lambda range")
    x = root.winfo_x()
    y = root.winfo_y()
    toplambda.geometry("+%d+%d" % (x + 500, y + 500))
    toplambda.focus_force()
    # toplambda.geometry("300x200")
    label1 = Label(toplambda, text="Start value")
    label2 = Label(toplambda, text="End value")
    label3 = Label(toplambda, text="Number of points")
    label1.grid(row=0, column=0)
    label2.grid(row=0, column=1)
    label3.grid(row=0, column=2)
    lambdaentry1 = Entry(toplambda)
    lambdaentry2 = Entry(toplambda)
    lambdaentry3 = Entry(toplambda)
    lambdaentry1.insert(0, "0.001")
    lambdaentry2.insert(0, "10")
    lambdaentry3.insert(0, "30")
    lambdaentry1.grid(row=1, column=0)
    lambdaentry2.grid(row=1, column=1)
    lambdaentry3.grid(row=1, column=2)

    def okcommand():
        global lamalist
        lamalist = np.geomspace(float(lambdaentry1.get()), float(lambdaentry2.get()), int(lambdaentry3.get()))
        toplambda.destroy()
        x = import_tree2.selection()
        datav_list = [datav_list_all[int(item_i)] for item_i in x]
        pointslist = DRTsolver(datav_list[0]).multilambda()
        amps = DRTsolver(datav_list[0]).amplitudes()
        # fig = plt.figure(3, figsize=(6, 4), dpi=80)
        fig, ax = plt.subplots()
        resnorm = [i[0] for i in pointslist]
        resnorm0 = pointslist[0][0]
        resnormrel = [i / resnorm0 for i in resnorm]
        solnorm = [i[1] for i in pointslist]
        solnnorm_n = [i / sum(amps) for i in solnorm]
        dat, = plt.plot(lamalist, solnnorm_n, 'o-', picker=8)
        annot = ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        plt.xscale("log")
        plt.ylim(0, 1)
        # plt.yscale("log")
        plt.xlabel("Lambda")
        plt.ylabel("Residual norm")
        names = [str(round(lamalist[n], 3)) for n in range(len(lamalist))]

        def update_annot(ind):
            x, y = dat.get_data()
            annot.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
            text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
            annot.set_text(text)
            annot.get_bbox_patch().set_alpha(0.4)

        def hover(event):
            vis = annot.get_visible()
            if event.inaxes == ax:
                cont, ind = dat.contains(event)
                if cont:
                    update_annot(ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

        def onpick1(event):
            if isinstance(event.artist, Line2D):
                thisline = event.artist
                xdata = thisline.get_xdata()
                ydata = thisline.get_ydata()
                ind = event.ind
                print('onpick1 line:', np.column_stack([xdata[ind], ydata[ind]]))

        # fig.canvas.mpl_connect('pick_event', onpick1)

        fig.canvas.mpl_connect("motion_notify_event", hover)

        plt.show()
        # except: messagebox.showerror("Error", "No data imported")

    okbutton = Button(toplambda, text="OK", command=okcommand)
    okbutton.grid(row=2, column=1)

    cancelbutton = Button(toplambda, text="Cancel", command=toplambda.destroy)
    cancelbutton.grid(row=2, column=2)

def normalize():
    plotselection2(0)

def lower_upperboundary():
    x = import_tree2.selection()
    datav_list = [datav_list_all[int(item_i)] for item_i in x]
    for i in x:
        instlist[int(i)].lowerhboundary()
    plotselection2(0)

def higher_upperboundary():
    x = import_tree2.selection()
    datav_list = [datav_list_all[int(item_i)] for item_i in x]
    for i in x:
        instlist[int(i)].higherhboundary()
    plotselection2(0)

def lower_lowerboundary():
    x = import_tree2.selection()
    datav_list = [datav_list_all[int(item_i)] for item_i in x]
    for i in x:
        instlist[int(i)].lowerlboundary()
    plotselection2(0)

def higher_lowerboundary():
    x = import_tree2.selection()
    datav_list = [datav_list_all[int(item_i)] for item_i in x]
    for i in x:
        instlist[int(i)].higherlboundary()
    plotselection2(0)

def getfminmax():
    if instlist:
        minmax = instlist[0].fminmax()

    else:
        minmax = (0,0)
    min_string = "{:.1E}".format(Decimal(str(round(minmax[0],1))))
    max_string = "{:.1E}".format(Decimal(str(round(minmax[1],1))))
    fmin_value.config(text=min_string+ " Hz")
    fmax_value.config(text=max_string+ " Hz")

def boundaryconfirm():
    x = import_tree2.selection()
    lvalue = float(fmin_value.get())
    hvalue = float(fmax_value.get())
    for i in x:
        instlist[int(i)].findclosest(lvalue, hvalue)

    plotselection2(0)

def draw_area():
    for i in range(50):
        pw.canvas.mpl_disconnect(i)
    if is_on:
        switch()

    draw_area_button.configure(bg="light green", fg="black")
    var_base.set(0)
    var_points.set(0)
    cid = pw.canvas.mpl_connect('button_press_event', pw.firstclick)
    cid2 = pw.canvas.mpl_connect('button_release_event', pw.secondclick)
    confirm_area_button.grid(row=1, column=1, padx=5, pady=0, sticky=W)
    #plotselection0(0)

def clearpoints():
    pw.clearallpoints()
    plotselection3(0)

def areas_calc(v1,v2,v3):
    values = slider_range.getValues()
    pw.plotrectangle(values[0], values[1])
    #print(values)

    #quickfit()

### DRT tab ###

imp2_frame = Frame(optitab, bg='grey')
imptree_frame = Frame(imp2_frame, bg='blue')
fig_frame2 = Frame(optitab)
settings_frame = Frame(optitab, bg='grey')
area_frame = Frame(settings_frame)
lama_frame = Frame(settings_frame)
t_constant_frame = Frame(settings_frame)
boundary_frame1 = Frame(settings_frame)
drt_right_frame = Frame(optitab, bg='grey')
drt_right_frame2 = Frame(optitab, bg='grey')
rcl_frame = Frame(fig_frame2, bg="grey")
legend_frame = Frame(fig_frame2, bg="grey")
boundary_frame = Frame(optitab, bg="grey")



# Define import tree 2
import_tree2 = ttk.Treeview(imptree_frame, selectmode='extended', show='tree')
vsbi2 = ttk.Scrollbar(imptree_frame, orient='vertical', command=import_tree2.yview)
vsbi2.pack(side='right', fill='y')
import_tree2.pack(side='left')
import_tree2.configure(yscrollcommand=vsbi2.set)
itemn2 = 0

l1 = Label(lama_frame, text="Choose lambda")
e1 = Entry(lama_frame)
settings_label = Label(settings_frame, text="Settings", font=('helvetica', 12, 'bold'))
exportentry = Entry(optitab)
e1.insert(0, 0.2)

arealabel = Label(area_frame, text="Area/cm²")
area_entry = Entry(area_frame)
area_entry.insert(0,1)
area_confirm = tk.Button(area_frame,text="Confirm", command=areaconfirm)

import_button2 = tk.Button(imp2_frame, text="Import data", command=multiimport2,font=('helvetica', 12, 'bold') )
import_settingsbutton2 = tk.Button(imp2_frame, text="Import settings", command=import_settings, font=('helvetica', 12, 'bold'))
exportbutton = tk.Button(drt_right_frame2, text="Export spectrum to .txt-file", command=savetxt, font=('helvetica', 12, 'bold'))
optibutton = tk.Button(drt_right_frame2, text="L-curve plot", command=lcurveplot, font=('helvetica', 12, 'bold'))
residualbutton = tk.Button(drt_right_frame2, text="Residual plot", command=residualplot, font=('helvetica', 12, 'bold'))

# Settings frame
# Number of time constants
options = [
    "1 N_f",
    "2 N_f",
    "3 N_f",
    "4 N_f",
    "6 N_f"
]
clicked = StringVar()
clicked.set(options[2])
drop = ttk.Combobox(t_constant_frame, value=options,state="readonly")
drop.current(2)
drop.bind("<<ComboboxSelected>>", plotselection2)
dropdownlabel = Label(t_constant_frame, text="Number of time constants")

boundaryoptions = ["+- 0 decade ", "+- 1 decade", "+- 2 decade", "+- 3 decade"]
clicked_b = StringVar()
clicked_b.set(boundaryoptions[1])
drop_b = ttk.Combobox(boundary_frame1, value=boundaryoptions,state="readonly")
drop_b.current(1)
drop_b.bind("<<ComboboxSelected>>", plotselection2)
dropdownlabel_b = Label(boundary_frame1, text="Boundaries")

var_r = IntVar()
var_l = IntVar()
var_c = IntVar()
var_r.set(1)
var_l.set(1)

gobutton = tk.Button(lama_frame, text="Go!", command=lambda: plotselection2(0))
checkframe = tk.Frame(settings_frame)
checklabel = Label(checkframe, text="Additional passive elements")
check_r = Checkbutton(checkframe, text="R", variable=var_r, command=lambda:plotselection2(0))
check_c = Checkbutton(checkframe, text="C", variable=var_c, command=lambda:plotselection2(0))
check_l = Checkbutton(checkframe, text="L", variable=var_l, command=lambda:plotselection2(0))


# DRT normalization
var_n = IntVar()
checkframe2 = tk.Frame(settings_frame)
check_n = Checkbutton(checkframe2, text="Relative Intensity", variable=var_n, command=normalize)

# RCL values
r_label = Label(rcl_frame, text="R")
c_label = Label(rcl_frame, text="C")
l_label = Label(rcl_frame, text="L")

# legend check
var_legend = IntVar()
check_legend = Checkbutton(legend_frame, text="Hide Legend", variable=var_legend, command=lambda: plotselection2(0))
check_legend.pack()

# upper lower boundaries
fmin_label = Label(boundary_frame, text="f_min", font=10)
fmin_value = Entry(boundary_frame)
fmax_label = Label(boundary_frame, text="f_max",font=10)
fmax_value = Entry(boundary_frame)

fmin_value.insert(0,0.01)
fmax_value.insert(0,1000000)

boundary_confirm_button = tk.Button(boundary_frame,text="Confirm", command=boundaryconfirm)


imp2_frame.grid(row=0, column=0, padx=5, pady=5,sticky=NW)
settings_frame.grid(row=1, column=0,rowspan=2, padx=5, pady=5,sticky=NW)
import_button2.grid(row=0, column=0, padx=10, pady=5)
import_settingsbutton2.grid(row=0, column=1, padx=10, pady=5)
imptree_frame.grid(row=2, column=0,columnspan=2, padx=10, pady=5,sticky=NW)

exportbutton.grid(row=0, column=0, padx=10, pady=10, sticky=NW)
optibutton.grid(row=1, column=0, padx=10, pady=10, sticky=NW)
residualbutton.grid(row=2, column=0, padx=10, pady=10, sticky=NW)
fig_frame2.grid(row=0, column=1, rowspan=3, padx=5, pady=5, sticky=NW)
drt_right_frame.grid(row=0, column=2, padx=5, pady=5, sticky=NW)
drt_right_frame2.grid(row=2, column=2, padx=5, pady=5, sticky=NW)
rcl_frame.pack(side=BOTTOM, pady=5)
legend_frame.pack(side=BOTTOM, pady=5)
boundary_frame.grid(row=1, column=2, padx=5, pady=5,sticky=NW)
# settings_frame
settings_label.grid(row=0, column=0, columnspan=2, padx=10, pady=5)
area_frame.grid(row=1, column=0,columnspan=2, padx=10, pady=10)
t_constant_frame.grid(row=2, column=0, columnspan=2, sticky=W, padx=10, pady=10)
boundary_frame1.grid(row=3, column=0, columnspan=2, sticky=W, padx=10, pady=10)
lama_frame.grid(row=5, column=0, columnspan=2, sticky=W, padx=10, pady=10)
checkframe.grid(row=4, column=0, sticky=W, padx=10, pady=10)
checkframe2.grid(row=4, column=1, padx=10, pady=10)

#normalize to area
var_area = IntVar()
check_area = Checkbutton(checkframe2, text="Relative intensity", variable=var_area, command=lambda:plotselection2(0))

create_canvas2(fig_frame2,(130,130))
create_canvas4(drt_right_frame, (90,60))

# tau frame
drop.grid(row=0, column=1, padx=10, pady=10)
dropdownlabel.grid(row=0, column=0, padx=10, pady=10)

# boundary frame
drop_b.grid(row=0, column=1, padx=10, pady=10)
dropdownlabel_b.grid(row=0, column=0, padx=10, pady=10)

# checkframe
checklabel.grid(row=0, column=0, columnspan=3)
check_r.grid(row=1, column=0, padx=10, pady=10)
check_c.grid(row=1, column=1, padx=10, pady=10)
check_l.grid(row=1, column=2, padx=10, pady=10)

# checkframe 2
check_n.grid(row=0, column=0, padx=10, pady=10)

# RCL label framel
r_label.grid(row=0, column=0, padx=10, pady=10)
c_label.grid(row=0, column=1, padx=10, pady=10)
l_label.grid(row=0, column=2, padx=10, pady=10)


# checkframe2
#check_n.grid(row=1, column=0, padx=10, pady=10)

# lama frame
l1.grid(row=0, column=0, padx=10, pady=10)
e1.grid(row=0, column=1, padx=10, pady=10)
gobutton.grid(row=0, column=2, pady=5, padx=5)

arealabel.grid(row=0, column=0, padx=10, pady=10)
area_entry.grid(row=0, column=1, padx=10, pady=10)
area_confirm.grid(row=0,column=2, padx=10, pady=10)

# boundary frame

#lower_upperboundary_button.grid(row=1, column=0, padx=10, pady=10)
#higher_upperboundary_button.grid(row=1, column=2, padx=10, pady=10)
#lower_lowerboundary_button.grid(row=1, column=3, padx=10, pady=10)
#higher_lowerboundary_button.grid(row=1, column=5, padx=10, pady=10)
fmin_value.grid(row=1, column=4, padx=10, pady=10)
fmax_value.grid(row=1, column=1, padx=10, pady=10)
fmin_label.grid(row=0,column=4, padx=5, pady=5, sticky=N)
fmax_label.grid(row=0,column=1, padx=5, pady=5, sticky=N)
boundary_confirm_button.grid(row=1, column=6, padx=10, pady=10)

### Fittab ###

fitframe_left = Frame(fittab, bg='grey')
drt_tree_frame = Frame(fitframe_left, bg='light blue')
sliderframe = Frame(fitframe_left)
fig_frame_fit = Frame(fittab, bg='grey')

var_a = IntVar()
var_a.set(1)
check_00 = Checkbutton(fitframe_left, text="Auto-detect peaks", variable=var_a, command=lambda: plotselection3(0))

var_d = IntVar()
var_d.set(1)
check_01 = Checkbutton(fitframe_left, text="Choose peaks manually", variable=var_d, command=lambda: plotselection3(0))
# Define import drt tree
import_drt_tree = ttk.Treeview(drt_tree_frame, selectmode='extended', show='tree')
vsbi4 = ttk.Scrollbar(drt_tree_frame, orient='vertical', command=import_drt_tree.yview)
vsbi4.pack(side='right', fill='y')
import_drt_tree.pack(side='left')
import_drt_tree.configure(yscrollcommand=vsbi4.set)
itemn4 = 0

pw = Plotwindow(fig_frame_fit,(130, 130))

import_drt_button = tk.Button(fitframe_left, text="Import DRT data", command=multiimport_drt,font=('helvetica', 12, 'bold') )
quickfitbutton = tk.Button(fitframe_left, text="Quickfit", command=quickfit,font=('helvetica', 12, 'bold') )
areafit_button = tk.Button(fitframe_left, text="Area fit", command=areafit,font=('helvetica', 12, 'bold') )
clearpoint_button = tk.Button(fitframe_left, text="Clear points", command=clearpoints,font=('helvetica', 12, 'bold') )
area_save_button = tk.Button(fitframe_left, text="Area save", command=lambda: savetxt_area(0),font=('helvetica', 12, 'bold') )

hleft = tk.DoubleVar(value=-6)
hright = tk.DoubleVar(value=1)
#slider_areas = Scale(sliderframe, orient=HORIZONTAL, resolution=0.01, from_=0.0, to=1.0, command=areas_calc)
slider_range = RangeSliderH(sliderframe, [hleft,hright], min_val=-10, max_val=4, padX=20, Width=200, Height=40, bar_radius=5, font_size=10)
slider_range.grid(row=0, column=0, padx=10, pady=10, sticky=NW)
hleft.trace_add('write',areas_calc)
hright.trace_add('write',areas_calc)


fitframe_left.grid(row=0, column=0, padx=10, pady=10, sticky=NW)
fig_frame_fit.grid(row=0, column=1, padx=10, pady=10, sticky=NW)
drt_tree_frame.grid(row=1, column=0, padx=5, pady=5, sticky=NW)
import_drt_button.grid(row=0,  column=0, padx=10, pady=10, sticky=NW)
quickfitbutton.grid(row=2,  column=0, padx=10, pady=10, sticky=NW)
clearpoint_button.grid(row=3,  column=0, padx=10, pady=10, sticky=NW)
sliderframe.grid(row=4,  column=0, padx=10, pady=10, sticky=NW)
#areafit_button.grid(row=2, column=0, padx=10, pady=10, sticky=NW)
#area_save_button.grid(row=4, column=0, padx=10, pady=10, sticky=NW)

# Mainloop
root.mainloop()
