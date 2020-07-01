import numpy as np

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

Flag = True

def _select_style(type, vc1, vc2, __warning, func1, func2, *arg, **kw):
    if vc1.get() == '固定' and vc2.get() == '固定':
        __warning('不被支持的原子选择方式')
        return
    '''
    if vc1.get() == '可调（仅原子）' and vc2.get() == '固定':
        func1(arg[-1], *arg[1:-1], arg[0], kw['raw']+kw['delta'])
    elif vc1.get() == '固定' and vc2.get() == '可调（仅原子）':
        func1(arg[0], *arg[1:-1], arg[-1], kw['raw']+kw['delta'])
    elif vc1.get() == '可调（仅原子）' and vc2.get() == '可调（仅原子）':
        func1(arg[-1], *arg[1:-1], arg[0], (kw['raw']+kw['delta']/2)*(-1)**type)
        func1(arg[0], *arg[1:-1], arg[-1], kw['raw']+kw['delta'])
    '''

    if vc1.get() == '可调（仅原子）' and vc2.get() == '固定':
        func1(arg[-1], *arg[1:-1], arg[0], x = kw['raw']+kw['delta'], mul = 1, dmul= 1, **kw)
    elif vc1.get() == '固定' and vc2.get() == '可调（仅原子）':
        func1(arg[0], *arg[1:-1], arg[-1], x = kw['raw']+kw['delta'], mul = 1, dmul=-1, **kw)
    elif vc1.get() == '可调（仅原子）' and vc2.get() == '可调（仅原子）':
        func1(arg[-1], *arg[1:-1], arg[0], x = kw['raw']+kw['delta']/2, mul = 0.5, dmul= 0.5, **kw)
        func1(arg[0], *arg[1:-1], arg[-1], x = kw['raw']+kw['delta'], mul = 1, dmul=-0.5, **kw)

    elif vc1.get() == '可调（基团）' and vc2.get() == '固定':
        func2(arg[-1], *arg[1:-1], arg[0], x = kw['raw']+kw['delta'], mul = 1, dmul= 1, **kw)
    elif vc1.get() == '固定' and vc2.get() == '可调（基团）':
        func2(arg[0], *arg[1:-1], arg[-1], x = kw['raw']+kw['delta'], mul = 1, dmul=-1, **kw)
    elif vc1.get() == '可调（基团）' and vc2.get() == '可调（基团）':
        func2(arg[-1], *arg[1:-1], arg[0], x = kw['raw']+kw['delta']/2, mul = 0.5, dmul= 0.5, **kw)
        func2(arg[0], *arg[1:-1], arg[-1], x = kw['raw']+kw['delta'], mul = 1, dmul=-0.5, **kw)
    else:
        __warning('暂未支持的原子选择方式')
        return


class _Bond_windows(tk.Toplevel):

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.wm_title("键设定")
        self.resizable(width=False, height=False)
        self.Molecule = self.parent.Molecule.copy()
        self.__select = []
        self.protocol("WM_DELETE_WINDOW", self.__quit)
        
        self.a = tk.IntVar()
        self.b = tk.IntVar()
        self.l = tk.DoubleVar()
        self.bl = tk.DoubleVar()
        self.l.trace("w", lambda name, index, mode: self.__change_text())
        self.bl.trace("w", lambda name, index, mode: self.__change_bl())
        self.vc1 = tk.StringVar()
        self.vc2 = tk.StringVar()
        style = ['固定', '可调（仅原子）', '可调（基团）']

        f0 = ttk.LabelFrame(self, padding=5, text='选择'); f0.grid(row=0, column=0)
        l1 = ttk.Label(f0, text='第一个原子的序号：'); l1.grid(row=0, column=0)
        l2 = ttk.Label(f0, text='第二个原子的序号：'); l2.grid(row=1, column=0)
        e1 = ttk.Entry(f0, textvariable=self.a); e1.grid(row=0, column=1)
        e2 = ttk.Entry(f0, textvariable=self.b); e2.grid(row=1, column=1)
        c1 = ttk.Combobox(f0, textvariable=self.vc1, values=style, exportselection=0, state='readonly')
        c1.grid(row=0, column=2); self.vc1.set('可调（仅原子）')
        c2 = ttk.Combobox(f0, textvariable=self.vc2, values=style, exportselection=0, state='readonly')
        c2.grid(row=1, column=2); self.vc2.set('可调（仅原子）')
        l4 = ttk.Label(f0, text='可在3D图形中点选'); l4.grid(row=2, columnspan=3)

        f1 = ttk.LabelFrame(self, padding=5, text='键级'); f1.grid(row=1, column=0)
        r1 = ttk.Radiobutton(f1, width=5, text='0', variable=self.bl, value=0.0); r1.grid(row=0, column=0)
        r2 = ttk.Radiobutton(f1, width=5, text='1', variable=self.bl, value=1.0); r2.grid(row=0, column=1)
        r3 = ttk.Radiobutton(f1, width=5, text='1.5', variable=self.bl, value=1.5); r3.grid(row=0, column=2)
        r4 = ttk.Radiobutton(f1, width=5, text='2', variable=self.bl, value=2.0); r4.grid(row=0, column=3)
        r5 = ttk.Radiobutton(f1, width=5, text='3', variable=self.bl, value=3.0); r5.grid(row=0, column=4)

        f2 = ttk.LabelFrame(self, padding=5, text='原子间距'); f2.grid(row=2, column=0)
        l3 = ttk.Label(f2, text='原子间距：'); l3.grid(row=0, column=0)
        e3 = ttk.Entry(f2, textvariable=self.l); e3.grid(row=0, column=1)
        self.s1 = tk.Scale(f2, resolution=0.001, from_=0.3, to=4, orient=tk.HORIZONTAL, length=250, command=self.__change_scale)
        self.s1.grid(row=1, columnspan=2)

        f3 = ttk.Frame(self, padding=5); f3.grid(row=3, column=0)
        b1 = ttk.Button(f3, text='确定', command=self.__commit); b1.grid(row=0, column=0)
        b2 = ttk.Button(f3, text='取消', command=self.__quit); b2.grid(row=0, column=1)

    def __modify_bond_length(self, Molecule):
        try:
            raw = Molecule.get_bond_length(self.a.get(), self.b.get())
            _select_style(0, self.vc1, self.vc2, self.__warning, Molecule.modify_bond_length_A, 
                        Molecule.modify_bond_length_G, self.a.get(), self.b.get(),
                        raw = raw, delta = self.l.get()-raw)
            self.parent.ddd.re_plot(self.Molecule)
        except RuntimeWarning:
            self.__warning('RuntimeWarning，请检查输入')

    def __change_bl(self):
        if self.bl.get() not in [0.0, 1.0, 1.5, 2.0, 3.0]: return
        self.Molecule.modify_bond_level(self.a.get(), self.b.get(), self.bl.get())
        self.parent.ddd.re_plot(self.Molecule)
        self.s1.set(self.l.get())

    def __change_scale(self, event):
        if '{:.3f}'.format(self.l.get()) != event:
            self.l.set(event)

    def __change_text(self):
        self.__modify_bond_length(self.Molecule)
        self.s1.set(self.l.get())

    def __commit(self):
        self.__modify_bond_length(self.parent.Molecule)
        self.__quit()
    
    def __warning(self, t):
        global Flag
        if Flag:
            Flag = False
            messagebox.showinfo('警告', t)   
            def ___a():
                global Flag
                Flag = True
            self.after(50, ___a())

    def __quit(self):
        self.parent.ddd.plot.clear_high_light()
        self.parent.ddd.re_plot(self.parent.Molecule)
        del self.parent._blw
        self.destroy()
    
    def __set_init(self):
        self.bl.set(self.parent.Molecule.get_bond_level(self.a.get(), self.b.get()))
        self.l.set(self.parent.Molecule.get_bond_length(self.a.get(), self.b.get()))
        self.s1.set(self.l.get())

    def select(self, i):
        self.__select.append(i)
        if len(self.__select) == 1:
            self.a.set(i)
        elif len(self.__select) == 2:
            self.b.set(i)
            self.__set_init()
        elif len(self.__select) >= 2:
            self.a.set(self.__select[-2])
            self.b.set(self.__select[-1])
            self.__set_init()
    
    def clear_selection(self):
        self.__select = []
        self.a.set(0)
        self.b.set(0)
        self.l.set(0)
        self.bl.set(0)

class _Bond_angle_windows(tk.Toplevel):

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.wm_title("键角")
        self.resizable(width=False, height=False)
        self.Molecule = self.parent.Molecule.copy()
        self.__select = []
        self.protocol("WM_DELETE_WINDOW", self.__quit)

        self.o = tk.IntVar()
        self.a = tk.IntVar()
        self.b = tk.IntVar()
        self.angle = tk.DoubleVar()
        self.angle.trace("w", lambda name, index, mode: self.__change_text())
        self.vc1 = tk.StringVar()
        self.vc3 = tk.StringVar()
        style = ['固定', '可调（仅原子）', '可调（基团）']

        f0 = ttk.LabelFrame(self, padding=5, text='选择'); f0.grid(row=0, column=0)
        l1 = ttk.Label(f0, text='第一个原子的序号：'); l1.grid(row=0, column=0)
        l2 = ttk.Label(f0, text='第二个原子的序号：'); l2.grid(row=1, column=0)
        l3 = ttk.Label(f0, text='第三个原子的序号：'); l3.grid(row=2, column=0)
        e1 = ttk.Entry(f0, textvariable=self.a); e1.grid(row=0, column=1)
        e2 = ttk.Entry(f0, textvariable=self.o); e2.grid(row=1, column=1)
        e3 = ttk.Entry(f0, textvariable=self.b); e3.grid(row=2, column=1)
        c1 = ttk.Combobox(f0, textvariable=self.vc1, values=style, exportselection=0, state='readonly')
        c1.grid(row=0, column=2); self.vc1.set('可调（仅原子）')
        c3 = ttk.Combobox(f0, textvariable=self.vc3, values=style, exportselection=0, state='readonly')
        c3.grid(row=2, column=2); self.vc3.set('可调（仅原子）')
        l5 = ttk.Label(f0, text='可在3D图形中点选'); l5.grid(row=3, columnspan=3)

        f2 = ttk.LabelFrame(self, padding=5, text='键角'); f2.grid(row=1, column=0)
        l4 = ttk.Label(f2, text='键角：'); l4.grid(row=0, column=0)
        e4 = ttk.Entry(f2, textvariable=self.angle); e4.grid(row=0, column=1)
        self.s1 = tk.Scale(f2, label='键角:', resolution=0.1, from_=1, to=180, orient=tk.HORIZONTAL, length=250, command=self.__change_scale)
        self.s1.grid(row=1, columnspan=2)

        f3 = ttk.Frame(self, padding=5); f3.grid(row=2, column=0)
        b1 = ttk.Button(f3, text='确定', command=self.__commit); b1.grid(row=4, column=0)
        b2 = ttk.Button(f3, text='取消', command=self.__quit); b2.grid(row=4, column=1)

    def __modify_bond_angle(self, Molecule):
        try:
            raw = Molecule.get_bond_angle(self.a.get(), self.o.get(), self.b.get(), 0)
            _select_style(0, self.vc1, self.vc3, self.__warning, Molecule.modify_bond_angle_A, 
                        Molecule.modify_bond_angle_G, self.a.get(), self.o.get(), self.b.get(), 
                        raw = raw, delta = self.angle.get()/180*np.pi - raw)
            self.parent.ddd.re_plot(self.Molecule)
        except RuntimeWarning:
            self.__warning('RuntimeWarning，请检查输入')

    def __change_scale(self, event):
        if '{:.1f}'.format(self.angle.get()) != event:
            self.angle.set(event)

    def __change_text(self):
        self.__modify_bond_angle(self.Molecule)
        self.s1.set(self.angle.get())

    def __commit(self):
        self.__modify_bond_angle(self.parent.Molecule)
        self.__quit()

    def __warning(self, t):
        global Flag
        if Flag:
            Flag = False
            messagebox.showinfo('警告', t)   
            def ___a():
                global Flag
                Flag = True
            self.after(50, ___a())
    
    def __quit(self):
        self.parent.ddd.plot.clear_high_light()
        self.parent.ddd.re_plot(self.parent.Molecule)
        del self.parent._baw
        self.destroy()

    def __set_init(self):
        a = self.parent.Molecule.get_bond_angle(self.a.get(), self.o.get(), self.b.get(), 0)/np.pi*180
        self.angle.set(a)
        self.s1.set(self.angle.get())

    def select(self, i):
        self.__select.append(i)
        if len(self.__select) == 1:
            self.a.set(i)
        elif len(self.__select) == 2:
            self.o.set(i)
        elif len(self.__select) == 3:
            self.b.set(i)
            self.__set_init()
        elif len(self.__select) >= 3:
            self.a.set(self.__select[-3])
            self.o.set(self.__select[-2])
            self.b.set(self.__select[-1])
            self.__set_init()
    
    def clear_selection(self):
        self.__select = []
        self.a.set(0)
        self.o.set(0)
        self.b.set(0)
        self.angle.set(0)


class _Dihedral_angle_windows(tk.Toplevel):

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.wm_title("二面角")
        self.resizable(width=False, height=False)
        self.Molecule = self.parent.Molecule.copy()
        self.__select = []
        self.protocol("WM_DELETE_WINDOW", self.__quit)

        self.a = tk.IntVar()
        self.b = tk.IntVar()
        self.c = tk.IntVar()
        self.d = tk.IntVar()
        self.angle = tk.DoubleVar()
        self.angle.trace("w", lambda name, index, mode: self.__change_text())
        self.vc1 = tk.StringVar()
        self.vc4 = tk.StringVar()
        style = ['固定', '可调（仅原子）', '可调（基团）']

        f0 = ttk.LabelFrame(self, padding=5, text='选择'); f0.grid(row=0, column=0)
        l1 = ttk.Label(f0, text='第一个原子的序号：'); l1.grid(row=0, column=0)
        l2 = ttk.Label(f0, text='第二个原子的序号：'); l2.grid(row=1, column=0)
        l3 = ttk.Label(f0, text='第三个原子的序号：'); l3.grid(row=2, column=0)
        l4 = ttk.Label(f0, text='第四个原子的序号：'); l4.grid(row=3, column=0)
        e1 = ttk.Entry(f0, textvariable=self.a); e1.grid(row=0, column=1)
        e2 = ttk.Entry(f0, textvariable=self.b); e2.grid(row=1, column=1)
        e3 = ttk.Entry(f0, textvariable=self.c); e3.grid(row=2, column=1)
        e4 = ttk.Entry(f0, textvariable=self.d); e4.grid(row=3, column=1)
        c1 = ttk.Combobox(f0, textvariable=self.vc1, values=style, exportselection=0, state='readonly')
        c1.grid(row=0, column=2); self.vc1.set('可调（仅原子）')
        c4 = ttk.Combobox(f0, textvariable=self.vc4, values=style, exportselection=0, state='readonly')
        c4.grid(row=3, column=2); self.vc4.set('可调（仅原子）')
        l6 = ttk.Label(f0, text='可在3D图形中点选'); l6.grid(row=4, columnspan=3)

        f2 = ttk.LabelFrame(self, padding=5, text='二面角'); f2.grid(row=1, column=0)
        l5 = ttk.Label(f2, text='二面角：'); l5.grid(row=0, column=0)
        e5 = ttk.Entry(f2, textvariable=self.angle); e5.grid(row=0, column=1)
        self.s1 = tk.Scale(f2, label='二面角:', resolution=0.1, from_=-180, to=180, orient=tk.HORIZONTAL, length=250, command=self.__change_scale)
        self.s1.grid(row=1, columnspan=2)

        f3 = ttk.Frame(self, padding=5); f3.grid(row=2, column=0)
        b1 = ttk.Button(f3, text='确定', command=self.__commit); b1.grid(row=0, column=0)
        b2 = ttk.Button(f3, text='取消', command=self.__quit); b2.grid(row=0, column=1)

    def __modify_dihedral_angle(self, Molecule):
        try:
            raw = Molecule.get_dihedral_angle(self.a.get(), self.b.get(), self.c.get(), self.d.get())
            _select_style(0, self.vc1, self.vc4, self.__warning, Molecule.modify_dihedral_angle_A, 
                        Molecule.modify_dihedral_angle_G, self.a.get(), self.b.get(), self.c.get(), self.d.get(), 
                        raw = raw, delta = self.angle.get()/180*np.pi - raw)
            self.parent.ddd.re_plot(self.Molecule)
        except RuntimeWarning:
            self.__warning('RuntimeWarning，请检查输入')

    def __change_scale(self, event):
        if '{:.1f}'.format(self.angle.get()) != event:
            self.angle.set(event)

    def __change_text(self):
        self.__modify_dihedral_angle(self.Molecule)
        self.s1.set(self.angle.get())

    def __commit(self):
        self.__modify_dihedral_angle(self.parent.Molecule)
        self.__quit()

    def __warning(self, t):
        global Flag
        if Flag:
            Flag = False
            messagebox.showinfo('警告', t)   
            def ___a():
                global Flag
                Flag = True
            self.after(50, ___a())

    def __quit(self):
        self.parent.ddd.plot.clear_high_light()
        self.parent.ddd.re_plot(self.parent.Molecule)
        del self.parent._daw
        self.destroy()

    def __set_init(self):
        self.angle.set(self.parent.Molecule.get_dihedral_angle(self.a.get(), self.b.get(), self.c.get(), self.d.get())/np.pi*180)
        self.s1.set(self.angle.get())

    def select(self, i):
        self.__select.append(i)
        if len(self.__select) == 1:
            self.a.set(i)
        elif len(self.__select) == 2:
            self.b.set(i)
        elif len(self.__select) == 3:
            self.c.set(i)
        elif len(self.__select) == 4:
            self.d.set(i)
            self.__set_init()
        elif len(self.__select) >= 4:
            self.a.set(self.__select[-4])
            self.b.set(self.__select[-3])
            self.c.set(self.__select[-2])
            self.d.set(self.__select[-1])
            self.__set_init()
    
    def clear_selection(self):
        self.__select = []
        self.a.set(0)
        self.b.set(0)
        self.c.set(0)
        self.d.set(0)
        self.angle.set(0)
