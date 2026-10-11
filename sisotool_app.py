import os
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import control as ct

plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

class EditLocationDialog(tk.Toplevel):
    def __init__(self, parent, title, initial_value):
        super().__init__(parent)
        self.title(title)
        self.geometry("440x160")
        self.resizable(False, False)
        self.transient(parent)
        self.result = None

        main_frame = ttk.Frame(self, padding=15)
        main_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        lbl = ttk.Label(main_frame, text="Enter new location (real or complex, e.g. -10 or -2+3j):")
        lbl.pack(side=tk.TOP, anchor="w", pady=(0, 8))
        
        self.ent = ttk.Entry(main_frame, width=45)
        self.ent.pack(side=tk.TOP, fill=tk.X, pady=(0, 15))
        self.ent.insert(0, str(initial_value))
        self.ent.select_range(0, tk.END)
        self.ent.focus()
        
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(side=tk.TOP, fill=tk.X)
        
        btn_ok = ttk.Button(btn_frame, text="OK", command=self.on_ok)
        btn_ok.pack(side=tk.RIGHT, padx=5)
        btn_cancel = ttk.Button(btn_frame, text="Cancel", command=self.destroy)
        btn_cancel.pack(side=tk.RIGHT, padx=5)
        
        self.ent.bind("<Return>", lambda e: self.on_ok())
        self.wait_visibility()
        try:
            self.grab_set()
        except Exception:
            pass

    def on_ok(self):
        self.result = self.ent.get()
        self.destroy()


class SISOToolApp:
    def __init__(self, root):
        self.root = root
        self.root.title("SISO Tool - Control System Designer (Python Edition)")
        self.root.geometry("1480x880")
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        
        self.G = ct.tf(0.01, [1, 12, 20])
        self.G_str = "0.01 / ((s+10)*(s+2))"
        
        self.C_zeros, self.C_poles, self.C_gain = [], [], 1.0
        self.C = ct.zpk(self.C_zeros, self.C_poles, self.C_gain)
        
        self.F_zeros, self.F_poles, self.F_gain = [], [], 1.0
        self.F = ct.zpk(self.F_zeros, self.F_poles, self.F_gain)
        
        self.C1_zeros, self.C1_poles, self.C1_gain = [], [], 1.0
        self.C1 = ct.zpk(self.C1_zeros, self.C1_poles, self.C1_gain)
        
        self.C2_zeros, self.C2_poles, self.C2_gain = [], [], 1.0
        self.C2 = ct.zpk(self.C2_zeros, self.C2_poles, self.C2_gain)
        
        self.H = ct.zpk([], [], 1)
        self.H_str = "1"
        self.G1 = ct.zpk([], [], 1)
        self.G1_str = "1"
        self.G2 = ct.zpk([], [], 1)
        self.G2_str = "1"
        self.Gd = ct.zpk([], [], 1)
        self.Gd_str = "1"
        self.H1 = ct.zpk([], [], 1)
        self.H1_str = "1"
        self.H2 = ct.zpk([], [], 1)
        self.H2_str = "1"
        
        self.architecture_idx = 1
        self.saved_designs = []
        self.saved_designs_meta = []
        self.flag_disp_C = []
        
        self.rl_mode = 'adjust'
        self.active_comp_name = "C"
        self.edit_sys = self.C
        self.custom_tfinal = ""

        print("[DEBUG] App Initialized with Plant G =", self.G_str)
        self.create_widgets()
        self.refresh_all()

    def on_close(self):
        print("[DEBUG] Application closing...")
        try:
            plt.close('all')
        except Exception:
            pass
        self.root.destroy()
        self.root.quit()

    def clean_val(self, val):
        if isinstance(val, (complex, np.complexfloating)):
            r = val.real
            i = val.imag
            if abs(r - round(r)) < 1e-4: r = round(r, 4)
            if abs(i - round(i)) < 1e-4: i = round(i, 4)
            return complex(r, i)
        else:
            v = float(val)
            if abs(v - round(v)) < 1e-4:
                v = round(v, 4)
            return v

    def format_loc(self, val):
        val = self.clean_val(val)
        if isinstance(val, (complex, np.complexfloating)):
            r, i = val.real, val.imag
            if abs(i) < 1e-5:
                return f"{r:g}"
            return f"{r:g}+{i:g}j" if i >= 0 else f"{r:g}{i:g}j"
        else:
            return f"{val:g}"

    def create_widgets(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.notebook.bind("<<NotebookTabChanged>>", self.on_tab_changed)
        
        self.tab_rl = ttk.Frame(self.notebook)
        self.tab_bode = ttk.Frame(self.notebook)
        self.tab_nyq = ttk.Frame(self.notebook)
        self.tab_all = ttk.Frame(self.notebook)
        self.tab_edit = ttk.Frame(self.notebook)
        
        self.notebook.add(self.tab_rl, text="Root Locus")
        self.notebook.add(self.tab_bode, text="Bode Diagram")
        self.notebook.add(self.tab_nyq, text="Nyquist")
        self.notebook.add(self.tab_all, text="All Plots Grid")
        self.notebook.add(self.tab_edit, text="Edit Controller")
        
        left_frame = ttk.Frame(self.root, width=520)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10, pady=10)
        
        self.fig_main, self.ax_main = plt.subplots(figsize=(5, 3.2))
        self.fig_main.patch.set_facecolor('white')
        self.canvas_main = FigureCanvasTkAgg(self.fig_main, master=left_frame)
        self.canvas_main.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        ctrl_panel = ttk.LabelFrame(left_frame, text="Design Controls")
        ctrl_panel.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=False, pady=5)
        
        ttk.Label(ctrl_panel, text="Plant G(s):").grid(row=0, column=0, sticky="w", padx=5, pady=4)
        self.ent_plant = ttk.Entry(ctrl_panel, width=28)
        self.ent_plant.insert(0, self.G_str)
        self.ent_plant.grid(row=0, column=1, columnspan=2, sticky="ew", padx=5, pady=4)
        self.ent_plant.bind("<Return>", self.update_plant)
        
        ttk.Label(ctrl_panel, text="View Response:").grid(row=1, column=0, sticky="w", padx=5, pady=4)
        self.cmb_response = ttk.Combobox(ctrl_panel, values=["Step Response", "Impulse Response", "Control Effort"], state="readonly")
        self.cmb_response.current(0)
        self.cmb_response.grid(row=1, column=1, columnspan=2, sticky="ew", padx=5, pady=4)
        self.cmb_response.bind("<<ComboboxSelected>>", lambda e: self.refresh_plots())
        
        ttk.Label(ctrl_panel, text="Sim Time (T_final):").grid(row=2, column=0, sticky="w", padx=5, pady=4)
        self.ent_tfinal = ttk.Entry(ctrl_panel, width=12)
        self.ent_tfinal.insert(0, "Auto")
        self.ent_tfinal.grid(row=2, column=1, sticky="ew", padx=5, pady=4)
        self.ent_tfinal.bind("<Return>", self.update_tfinal_from_box)
        
        self.btn_save = ttk.Button(ctrl_panel, text="Save", command=self.save_controller)
        self.btn_save.grid(row=3, column=0, sticky="ew", padx=2, pady=4)
        
        self.btn_edit_ctrl = ttk.Button(ctrl_panel, text="Edit", command=lambda: self.notebook.select(self.tab_edit))
        self.btn_edit_ctrl.grid(row=3, column=1, sticky="ew", padx=2, pady=4)
        
        self.btn_arch = ttk.Button(ctrl_panel, text="Architecture...", command=self.open_architecture_window)
        self.btn_arch.grid(row=3, column=2, sticky="ew", padx=2, pady=4)

        ttk.Label(ctrl_panel, text="Saved Designs (Double-click to Toggle):").grid(row=4, column=0, columnspan=3, sticky="w", padx=5, pady=2)
        
        list_frame = ttk.Frame(ctrl_panel)
        list_frame.grid(row=5, column=0, columnspan=3, sticky="ew", padx=5, pady=2)
        
        self.listbox_designs = tk.Listbox(list_frame, height=3, exportselection=False)
        self.listbox_designs.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        # Use Double-click to toggle visibility so single-click can select for actions
        self.listbox_designs.bind("<Double-1>", self.toggle_design_visibility)
        
        # Side-by-side action buttons frame
        btn_action_frame = ttk.Frame(ctrl_panel)
        btn_action_frame.grid(row=6, column=0, columnspan=3, sticky="ew", padx=2, pady=2)
        
        self.btn_make_active = ttk.Button(btn_action_frame, text="Make Active", command=self.make_active_design)
        self.btn_make_active.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2)
        
        self.btn_del_design = ttk.Button(btn_action_frame, text="Delete Selected", command=self.delete_selected_design)
        self.btn_del_design.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2)
        
        self.lbl_status = ttk.Label(ctrl_panel, text="Architecture: Config 1 | Status: Ready", font=("Arial", 9, "bold"), foreground="blue")
        self.lbl_status.grid(row=7, column=0, columnspan=3, sticky="w", padx=5, pady=2)
        
        ttk.Label(ctrl_panel, text="Gain (k):").grid(row=8, column=0, sticky="w", padx=5, pady=4)
        self.slider_gain = ttk.Scale(ctrl_panel, from_=0.01, to=2000, value=1, command=self.update_gain_slider)
        self.slider_gain.grid(row=8, column=1, sticky="ew", padx=5, pady=4)
        
        self.lbl_gain_val = ttk.Label(ctrl_panel, text="1.00", font=("Arial", 9, "bold"))
        self.lbl_gain_val.grid(row=8, column=2, sticky="w", padx=5, pady=4)
        
        self.setup_tab_plots()
        self.setup_edit_controller_ui()

    def update_tfinal_from_box(self, event=None):
        val = self.ent_tfinal.get().strip()
        if val.lower() in ["", "auto", "none"]:
            self.custom_tfinal = ""
            self.ent_tfinal.delete(0, tk.END)
            self.ent_tfinal.insert(0, "Auto")
        else:
            try:
                self.custom_tfinal = float(val)
            except ValueError:
                messagebox.showerror("Error", "Please enter a valid number for T_final or 'Auto'")
                return
        self.refresh_plots()

    def setup_tab_plots(self):
        rl_toolbar = ttk.Frame(self.tab_rl)
        rl_toolbar.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        
        ttk.Button(rl_toolbar, text="+ Real Pole", command=lambda: self.set_rl_mode('add_rpole')).pack(side=tk.LEFT, padx=2)
        ttk.Button(rl_toolbar, text="+ Complex Pole", command=lambda: self.set_rl_mode('add_cpole')).pack(side=tk.LEFT, padx=2)
        ttk.Button(rl_toolbar, text="+ Real Zero", command=lambda: self.set_rl_mode('add_rzero')).pack(side=tk.LEFT, padx=2)
        ttk.Button(rl_toolbar, text="+ Complex Zero", command=lambda: self.set_rl_mode('add_czero')).pack(side=tk.LEFT, padx=2)
        ttk.Button(rl_toolbar, text="Delete P/Z", command=lambda: self.set_rl_mode('delete')).pack(side=tk.LEFT, padx=2)
        ttk.Button(rl_toolbar, text="Move / Adjust", command=lambda: self.set_rl_mode('adjust')).pack(side=tk.LEFT, padx=2)
        
        self.lbl_mode_status = ttk.Label(rl_toolbar, text="Mode: Adjust / Drag", font=("Arial", 9, "italic"), foreground="green")
        self.lbl_mode_status.pack(side=tk.RIGHT, padx=5)

        self.fig_rl, self.ax_rl = plt.subplots(figsize=(6, 6))
        self.fig_rl.patch.set_facecolor('white')
        self.canvas_rl = FigureCanvasTkAgg(self.fig_rl, master=self.tab_rl)
        self.canvas_rl.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.fig_rl.canvas.mpl_connect('button_press_event', self.on_root_locus_click)
        
        self.fig_bode, (self.ax_bode_m, self.ax_bode_p) = plt.subplots(2, 1, figsize=(6, 6))
        self.fig_bode.patch.set_facecolor('white')
        self.canvas_bode = FigureCanvasTkAgg(self.fig_bode, master=self.tab_bode)
        self.canvas_bode.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.fig_nyq, self.ax_nyq = plt.subplots(figsize=(6, 6))
        self.fig_nyq.patch.set_facecolor('white')
        self.canvas_nyq = FigureCanvasTkAgg(self.fig_nyq, master=self.tab_nyq)
        self.canvas_nyq.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.fig_all, axs = plt.subplots(2, 2, figsize=(7, 7))
        self.fig_all.patch.set_facecolor('white')
        self.ax_all_rl, self.ax_all_bm = axs[0, 0], axs[0, 1]
        self.ax_all_ny, self.ax_all_bp = axs[1, 0], axs[1, 1]
        self.canvas_all = FigureCanvasTkAgg(self.fig_all, master=self.tab_all)
        self.canvas_all.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def setup_edit_controller_ui(self):
        top_edit_frame = ttk.LabelFrame(self.tab_edit, text="Compensator")
        top_edit_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)
        
        self.cmb_comp_select = ttk.Combobox(top_edit_frame, values=["C", "F"], state="readonly", width=5)
        self.cmb_comp_select.current(0)
        self.cmb_comp_select.grid(row=0, column=0, padx=8, pady=10, sticky="w")
        self.cmb_comp_select.bind("<<ComboboxSelected>>", self.on_compensator_selected)
        
        ttk.Label(top_edit_frame, text=" = ").grid(row=0, column=1)
        
        self.ent_comp_gain = ttk.Entry(top_edit_frame, width=12)
        self.ent_comp_gain.insert(0, "1.0")
        self.ent_comp_gain.grid(row=0, column=2, padx=5, pady=10, sticky="ew")
        self.ent_comp_gain.bind("<Return>", self.update_comp_gain_from_box)
        
        ttk.Label(top_edit_frame, text=" x ").grid(row=0, column=3)
        
        self.lbl_tf_display = ttk.Label(top_edit_frame, text="1", font=("Courier", 10), justify=tk.LEFT)
        self.lbl_tf_display.grid(row=0, column=4, padx=15, pady=10, sticky="w")
        
        table_frame = ttk.LabelFrame(self.tab_edit, text="Compensator Dynamics (Double-click Row to Edit Location)")
        table_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        columns = ("type", "location", "damping", "frequency")
        self.tree_pz = ttk.Treeview(table_frame, columns=columns, show="headings", height=10)
        self.tree_pz.heading("type", text="Type")
        self.tree_pz.heading("location", text="Location")
        self.tree_pz.heading("damping", text="Damping")
        self.tree_pz.heading("frequency", text="Frequency")
        self.tree_pz.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.tree_pz.bind("<Double-1>", self.on_table_row_double_click)
        
        scrollbar = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.tree_pz.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree_pz.configure(yscrollcommand=scrollbar.set)

    def set_rl_mode(self, mode):
        self.rl_mode = mode
        mode_names = {
            'adjust': 'Move / Adjust',
            'add_rpole': 'Add Real Pole',
            'add_cpole': 'Add Complex Pole',
            'add_rzero': 'Add Real Zero',
            'add_czero': 'Add Complex Zero',
            'delete': 'Delete Pole/Zero'
        }
        self.lbl_mode_status.config(text=f"Mode: {mode_names.get(mode, mode)}")

    def on_table_row_double_click(self, event):
        row_id = self.tree_pz.identify_row(event.y)
        if not row_id:
            return
        parts = row_id.split('_')
        if len(parts) != 2:
            return
        item_type, target_list_idx_str = parts[0], parts[1]
        target_list_idx = int(target_list_idx_str)
        is_pole = (item_type == "pole")

        if self.active_comp_name == "C":
            zeros, poles, k = self.C_zeros, self.C_poles, self.C_gain
        elif self.active_comp_name == "F":
            zeros, poles, k = self.F_zeros, self.F_poles, self.F_gain
        elif self.active_comp_name == "C1":
            zeros, poles, k = self.C1_zeros, self.C1_poles, self.C1_gain
        elif self.active_comp_name == "C2":
            zeros, poles, k = self.C2_zeros, self.C2_poles, self.C2_gain
        else:
            zeros, poles, k = self.C_zeros, self.C_poles, self.C_gain

        current_val = poles[target_list_idx] if is_pole else zeros[target_list_idx]
        dlg = EditLocationDialog(self.root, "Edit Location", current_val)
        self.root.wait_window(dlg)
        new_val_str = dlg.result

        if new_val_str is not None:
            try:
                expr = new_val_str.strip().replace('i', 'j')
                new_val = complex(expr) if 'j' in expr else float(expr)
                new_val = self.clean_val(new_val)

                if is_pole:
                    poles[target_list_idx] = new_val
                    if abs(np.imag(current_val)) > 1e-5 and abs(np.imag(new_val)) > 1e-5:
                        conj_target = np.conj(current_val)
                        for i, p in enumerate(poles):
                            if abs(p - conj_target) < 1e-4:
                                poles[i] = np.conj(new_val)
                                break
                else:
                    zeros[target_list_idx] = new_val
                    if abs(np.imag(current_val)) > 1e-5 and abs(np.imag(new_val)) > 1e-5:
                        conj_target = np.conj(current_val)
                        for i, z in enumerate(zeros):
                            if abs(z - conj_target) < 1e-4:
                                zeros[i] = np.conj(new_val)
                                break

                new_sys = ct.zpk(zeros, poles, k)
                if self.active_comp_name == "C":
                    self.C = new_sys
                    self.edit_sys = self.C
                elif self.active_comp_name == "F":
                    self.F = new_sys
                    self.edit_sys = self.F
                elif self.active_comp_name == "C1":
                    self.C1 = new_sys
                    self.edit_sys = self.C1
                elif self.active_comp_name == "C2":
                    self.C2 = new_sys
                    self.edit_sys = self.C2

                self.refresh_all()
            except Exception as e:
                messagebox.showerror("Error", f"Invalid numeric or complex value: {e}", parent=self.root)

    def on_root_locus_click(self, event):
        if event.inaxes != self.ax_rl:
            return
        x, y = event.xdata, event.ydata
        if x is None or y is None:
            return
            
        if self.active_comp_name == "C":
            zeros, poles, k = self.C_zeros, self.C_poles, self.C_gain
        elif self.active_comp_name == "F":
            zeros, poles, k = self.F_zeros, self.F_poles, self.F_gain
        elif self.active_comp_name == "C1":
            zeros, poles, k = self.C1_zeros, self.C1_poles, self.C1_gain
        elif self.active_comp_name == "C2":
            zeros, poles, k = self.C2_zeros, self.C2_poles, self.C2_gain
        else:
            zeros, poles, k = self.C_zeros, self.C_poles, self.C_gain

        if self.rl_mode == 'add_rpole':
            poles.append(self.clean_val(float(x)))
        elif self.rl_mode == 'add_cpole':
            poles.extend([self.clean_val(complex(x, abs(y))), self.clean_val(complex(x, -abs(y)))])
        elif self.rl_mode == 'add_rzero':
            zeros.append(self.clean_val(float(x)))
        elif self.rl_mode == 'add_czero':
            zeros.extend([self.clean_val(complex(x, abs(y))), self.clean_val(complex(x, -abs(y)))])
        elif self.rl_mode == 'delete':
            all_pts = [(p, 'pole', i) for i, p in enumerate(poles)] + [(z, 'zero', i) for i, z in enumerate(zeros)]
            if all_pts:
                distances = [np.hypot(pt[0].real - x, pt[0].imag - y) for pt in all_pts]
                min_idx = np.argmin(distances)
                pt, ptype, idx = all_pts[min_idx]
                if ptype == 'pole':
                    poles.pop(idx)
                    if abs(np.imag(pt)) > 1e-5:
                        conj_target = np.conj(pt)
                        for i, p in enumerate(poles):
                            if abs(p - conj_target) < 1e-4:
                                poles.pop(i)
                                break
                else:
                    zeros.pop(idx)
                    if abs(np.imag(pt)) > 1e-5:
                        conj_target = np.conj(pt)
                        for i, z in enumerate(zeros):
                            if abs(z - conj_target) < 1e-4:
                                zeros.pop(i)
                                break
        
        new_sys = ct.zpk(zeros, poles, k)
        if self.active_comp_name == "C":
            self.C = new_sys
            self.edit_sys = self.C
        elif self.active_comp_name == "F":
            self.F = new_sys
            self.edit_sys = self.F
        elif self.active_comp_name == "C1":
            self.C1 = new_sys
            self.edit_sys = self.C1
        elif self.active_comp_name == "C2":
            self.C2 = new_sys
            self.edit_sys = self.C2
            
        self.refresh_all()

    def safe_tf(self, sys_obj):
        """Safely convert any python object/number to an LTI transfer function object."""
        if isinstance(sys_obj, ct.LTI):
            return sys_obj
        try:
            if isinstance(sys_obj, (int, float, complex)):
                return ct.tf([sys_obj], [1])
            return ct.tf(sys_obj)
        except Exception:
            return ct.tf([1], [1])

    def get_active_system_eval(self, C_override=None):
        cfg = self.architecture_idx
        F = self.safe_tf(self.F)
        C = self.safe_tf(C_override if C_override is not None else self.C)
        C1 = self.safe_tf(self.C1)
        C2 = self.safe_tf(self.C2)
        G = self.safe_tf(self.G)
        G1 = self.safe_tf(self.G1)
        G2 = self.safe_tf(self.G2)
        H = self.safe_tf(self.H)
        
        if cfg == 1:
            sys_ol = C * G * H
            sys_cl = ct.feedback(F * C * G, H)
        elif cfg == 2:
            sys_ol = G * C * H
            sys_cl = ct.feedback(F * G, C * H)
        elif cfg == 3:
            sys_ol = C * G * H
            sys_cl = ct.feedback(F * G, C * H)
        elif cfg == 4:
            sys_ol = C1 * C2 * G * H
            sys_cl = ct.feedback(F * C1 * C2 * G, H)
        elif cfg == 5:
            sys_ol = C * (G1 + G2)
            sys_cl = ct.feedback(F * C * (G1 + G2), ct.tf(1, 1))
        elif cfg == 6:
            H1 = self.safe_tf(self.H1)
            sys_ol = C1 * G1 * H1
            sys_cl = ct.feedback(F * C1 * G1, H1)
        else:
            sys_ol = C * G * H
            sys_cl = ct.feedback(F * C * G, H)
            
        return sys_ol, sys_cl

    def update_plant(self, event=None):
        expr = self.ent_plant.get()
        try:
            s = ct.tf('s')
            self.G = eval(expr)
            self.G_str = expr
            self.lbl_status.config(text=f"Architecture: Config {self.architecture_idx} | Plant updated successfully.")
            self.refresh_all()
        except Exception as e:
            messagebox.showerror("Syntax Error", f"Invalid Transfer Function:\n{e}")

    def update_gain_slider(self, val):
        k = float(val)
        self.C_gain = k
        self.C = ct.zpk(self.C_zeros, self.C_poles, self.C_gain)
        self.lbl_gain_val.config(text=f"{k:.2f}")
        if self.active_comp_name == "C":
            self.edit_sys = self.C
        self.refresh_all()

    def update_comp_gain_from_box(self, event=None):
        try:
            k = float(self.ent_comp_gain.get())
            if self.active_comp_name == "C":
                self.C_gain = k
                self.C = ct.zpk(self.C_zeros, self.C_poles, self.C_gain)
                self.edit_sys = self.C
                if k > self.slider_gain.cget("to"):
                    self.slider_gain.configure(to=max(2000, k * 2))
                self.slider_gain.set(k)
                self.lbl_gain_val.config(text=f"{k:.2f}")
            elif self.active_comp_name == "F":
                self.F_gain = k
                self.F = ct.zpk(self.F_zeros, self.F_poles, self.F_gain)
                self.edit_sys = self.F
            elif self.active_comp_name == "C1":
                self.C1_gain = k
                self.C1 = ct.zpk(self.C1_zeros, self.C1_poles, self.C1_gain)
                self.edit_sys = self.C1
            elif self.active_comp_name == "C2":
                self.C2_gain = k
                self.C2 = ct.zpk(self.C2_zeros, self.C2_poles, self.C2_gain)
                self.edit_sys = self.C2
                
            self.refresh_all()
        except Exception as e:
            messagebox.showerror("Error", f"Invalid gain value: {e}")

    def on_compensator_selected(self, event=None):
        sel = self.cmb_comp_select.get()
        self.active_comp_name = sel
        if sel == "C":
            self.edit_sys = self.C
        elif sel == "F":
            self.edit_sys = self.F
        elif sel == "C1":
            self.edit_sys = self.C1
        elif sel == "C2":
            self.edit_sys = self.C2
        self.refresh_edit_controller_panel()

    def save_controller(self):
        self.saved_designs.append(self.C)
        self.saved_designs_meta.append(f"Design {len(self.saved_designs)}")
        self.flag_disp_C.append(1)
        self.update_design_listbox()
        messagebox.showinfo("Saved", f"Controller design saved! Total saved: {len(self.saved_designs)}")
        self.refresh_plots()

    def update_design_listbox(self):
        self.listbox_designs.delete(0, tk.END)
        for idx, name in enumerate(self.saved_designs_meta):
            status = "[Active]" if self.flag_disp_C[idx] == 1 else "[Hidden]"
            self.listbox_designs.insert(tk.END, f"{name} {status}")

    def toggle_design_visibility(self, event):
        sel = self.listbox_designs.curselection()
        if sel:
            idx = sel[0]
            self.flag_disp_C[idx] = 0 if self.flag_disp_C[idx] == 1 else 1
            self.update_design_listbox()
            self.refresh_plots()

    def delete_selected_design(self):
        sel = self.listbox_designs.curselection()
        if sel:
            idx = sel[0]
            del self.saved_designs[idx]
            del self.saved_designs_meta[idx]
            del self.flag_disp_C[idx]
            self.update_design_listbox()
            self.refresh_plots()
        else:
            messagebox.showwarning("Warning", "Please select a saved design from the list to delete.")

    def get_transfer_function_string(self, sys):
        try:
            num, den = ct.tfdata(sys)
            num_arr = np.array(num[0][0], dtype=float)
            den_arr = np.array(den[0][0], dtype=float)
            
            num_arr[abs(num_arr - np.round(num_arr)) < 1e-4] = np.round(num_arr[abs(num_arr - np.round(num_arr)) < 1e-4])
            den_arr[abs(den_arr - np.round(den_arr)) < 1e-4] = np.round(den_arr[abs(den_arr - np.round(den_arr)) < 1e-4])
            
            num_list = [float(x) for x in num_arr]
            den_list = [float(x) for x in den_arr]
            
            return f"ct.tf({num_list}, {den_list})"
        except Exception:
            return "1"

    def get_sys_string(self, name):
        if name in ['G', 'G1', 'G2', 'Gd', 'H1', 'H2']:
            return getattr(self, f"{name}_str", "1")
        if name in ['C', 'F', 'C1', 'C2']:
            sys_obj = getattr(self, name, None)
            if sys_obj is not None:
                return self.get_transfer_function_string(sys_obj)
        return "1"

    def open_architecture_window(self):
        win = ToplevelArchitecture(self)

    def refresh_edit_controller_panel(self):
        self.cmb_comp_select.set(self.active_comp_name)
        
        if self.active_comp_name == "C":
            k = self.C_gain
            poles = self.C_poles
            zeros = self.C_zeros
        elif self.active_comp_name == "F":
            k = self.F_gain
            poles = self.F_poles
            zeros = self.F_zeros
        elif self.active_comp_name == "C1":
            k = self.C1_gain
            poles = self.C1_poles
            zeros = self.C1_zeros
        elif self.active_comp_name == "C2":
            k = self.C2_gain
            poles = self.C2_poles
            zeros = self.C2_zeros
        else:
            k = 1.0
            poles, zeros = [], []
            
        if self.root.focus_get() != self.ent_comp_gain:
            self.ent_comp_gain.delete(0, tk.END)
            self.ent_comp_gain.insert(0, f"{k}")
        
        for row in self.tree_pz.get_children():
            self.tree_pz.delete(row)
            
        for i, p in enumerate(poles):
            p_clean = self.clean_val(p)
            damp = -np.cos(np.angle(p_clean)) if np.imag(p_clean) != 0 else 1.0
            freq = np.abs(p_clean)
            loc_str = self.format_loc(p_clean)
            type_str = "Complex Pole" if abs(np.imag(p_clean)) > 1e-5 else "Pole"
            self.tree_pz.insert("", tk.END, iid=f"pole_{i}", values=(type_str, loc_str, f"{damp:.3f}", f"{freq:.3f}"))
            
        for i, z in enumerate(zeros):
            z_clean = self.clean_val(z)
            damp = -np.cos(np.angle(z_clean)) if np.imag(z_clean) != 0 else 1.0
            freq = np.abs(z_clean)
            loc_str = self.format_loc(z_clean)
            type_str = "Complex Zero" if abs(np.imag(z_clean)) > 1e-5 else "Zero"
            self.tree_pz.insert("", tk.END, iid=f"zero_{i}", values=(type_str, loc_str, f"{damp:.3f}", f"{freq:.3f}"))
            
        try:
            zeros_list = list(self.edit_sys.zeros())
            poles_list = list(self.edit_sys.poles())
            unit_sys = ct.zpk(zeros_list, poles_list, 1.0)
            
            num_poly = unit_sys.num[0][0]
            den_poly = unit_sys.den[0][0]
            
            def poly_to_str(poly):
                n = len(poly)
                terms = []
                for idx, val in enumerate(poly):
                    power = n - 1 - idx
                    val_c = self.clean_val(val)
                    if abs(val_c) < 1e-12:
                        continue
                    v_abs = abs(val_c)
                    
                    if power > 1:
                        coeff_str = "" if abs(v_abs - 1) < 1e-12 else f"{v_abs:g}"
                        term = f"{coeff_str}s^{power}"
                    elif power == 1:
                        coeff_str = "" if abs(v_abs - 1) < 1e-12 else f"{v_abs:g}"
                        term = f"{coeff_str}s"
                    else:
                        term = f"{v_abs:g}"
                    terms.append((val_c, term))
                
                if not terms:
                    return "0"
                
                s = ""
                for i, (val_c, term) in enumerate(terms):
                    if i == 0:
                        sign = "-" if val_c < 0 else ""
                        s += sign + term
                    else:
                        sign = " - " if val_c < 0 else " + "
                        s += sign + term
                return s

            num_str = poly_to_str(num_poly)
            den_str = poly_to_str(den_poly)
            max_len = max(len(num_str), len(den_str))
            dash_str = "-" * max_len
            
            num_str = num_str.center(max_len)
            den_str = den_str.center(max_len)
            
            tf_clean_repr = f"{num_str}\n{dash_str}\n{den_str}"
            self.lbl_tf_display.config(text=tf_clean_repr)
        except Exception:
            self.lbl_tf_display.config(text=str(self.edit_sys))

    def on_tab_changed(self, event):
        if self.notebook.select() == str(self.tab_edit):
            self.refresh_edit_controller_panel()

    def refresh_all(self):
        self.refresh_plots()
        if self.notebook.select() == str(self.tab_edit):
            self.refresh_edit_controller_panel()

    def plot_bode_complete(self, sys, ax_m, ax_p, title_str, is_grid=False):
        ax_m.clear()
        ax_p.clear()
        
        omega = np.logspace(-1, 3, 500)
        fresp = ct.frequency_response(sys, omega)
        mag = fresp.magnitude
        phase = np.rad2deg(np.unwrap(fresp.phase))
        omega_pts = fresp.omega
        
        ax_m.semilogx(omega_pts, 20 * np.log10(mag))
        ax_m.set_title(title_str, fontsize=9 if is_grid else 10)
        ax_m.set_ylabel("Magnitude (dB)", fontsize=8 if is_grid else 10)
        ax_m.grid(True, which="both", ls="--")
        ax_m.set_xlim(omega_pts[0], omega_pts[-1])
        
        ax_p.semilogx(omega_pts, phase)
        ax_p.set_xlabel("Frequency (rad/s)", fontsize=8 if is_grid else 10)
        ax_p.set_ylabel("Phase (deg)", fontsize=8 if is_grid else 10)
        ax_p.grid(True, which="both", ls="--")
        ax_p.set_xlim(omega_pts[0], omega_pts[-1])
        
        try:
            gm, pm, sm, wg, wp, ws = ct.stability_margins(sys)
            gm_db = 20 * np.log10(gm) if gm > 0 and not np.isinf(gm) else float('inf')
            gm_str = f"G.M: {gm_db:.1f} dB" if not np.isinf(gm_db) else "G.M: Inf"
            pm_str = f"P.M: {pm:.1f} deg"
            wg_val = f"{wg:.2f}" if not np.isnan(wg) else "NaN"
            wp_val = f"{wp:.2f}" if not np.isnan(wp) else "NaN"
        except Exception:
            gm_str, pm_str, wg_val, wp_val = "G.M: N/A", "P.M: N/A", "NaN", "NaN"

        font_sz = 7 if is_grid else 8
        ax_m.text(0.05, 0.12, f"{gm_str}\nFreq: {wg_val}", transform=ax_m.transAxes, fontsize=font_sz, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax_p.text(0.05, 0.12, f"{pm_str}\nFreq: {wp_val}", transform=ax_p.transAxes, fontsize=font_sz, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        try:
            p = sys.poles()
            z = sys.zeros()
            
            if len(z) > 0:
                for z_val in z:
                    z_mag = abs(z_val)
                    idx = 0 if z_mag < 1e-5 else np.argmin(np.abs(omega_pts - z_mag))
                    ax_m.semilogx(omega_pts[idx], 20 * np.log10(mag[idx]), 'go', markersize=7, markeredgewidth=1.5, zorder=4)
                    ax_p.semilogx(omega_pts[idx], phase[idx], 'go', markersize=7, markeredgewidth=1.5, zorder=4)
            
            if len(p) > 0:
                for p_val in p:
                    p_mag = abs(p_val)
                    idx = 0 if p_mag < 1e-5 else np.argmin(np.abs(omega_pts - p_mag))
                    ax_m.semilogx(omega_pts[idx], 20 * np.log10(mag[idx]), 'rx', markersize=9, markeredgewidth=2.5, zorder=5)
                    ax_p.semilogx(omega_pts[idx], phase[idx], 'rx', markersize=9, markeredgewidth=2.5, zorder=5)
        except Exception:
            pass

    def draw_octave_style_nyquist(self, sys, ax, title_str):
        ax.clear()
        omega = np.logspace(-3, 3, 2000)
        fresp = ct.frequency_response(sys, omega)
        
        H = fresp.complex
        real_pos = np.real(H)
        imag_pos = np.imag(H)
        
        real_neg = real_pos[::-1]
        imag_neg = -imag_pos[::-1]
        
        real_full = np.concatenate([real_neg, real_pos])
        imag_full = np.concatenate([imag_neg, imag_pos])
        
        ax.plot(real_full, imag_full, 'b-', linewidth=1.5)
        ax.plot(-1, 0, 'r+', markersize=8, markeredgewidth=1.5)
        ax.set_title(title_str, fontsize=10)
        ax.set_xlabel("Real axis", fontsize=9)
        ax.set_ylabel("Imaginary axis", fontsize=9)
        ax.grid(True, linestyle='--', alpha=0.7)

    def refresh_plots(self):
        sys_ol, sys_cl = self.get_active_system_eval(self.C)

        t_vec = None
        if isinstance(self.custom_tfinal, (int, float)) and self.custom_tfinal > 0:
            t_vec = np.linspace(0, self.custom_tfinal, 500)

        try:
            self.ax_main.clear()
            choice = self.cmb_response.current()
            
            # Plot Active Design with a solid, prominent line
            if choice == 0:
                t, y = ct.step_response(sys_cl, T=t_vec)
                self.ax_main.plot(t, y, linewidth=2, label="Active (Current)")
                self.ax_main.set_title("Closed-Loop Step Response", fontsize=10)
            elif choice == 1:
                t, y = ct.impulse_response(sys_cl, T=t_vec)
                self.ax_main.plot(t, y, linewidth=2, label="Active (Current)")
                self.ax_main.set_title("Closed-Loop Impulse Response", fontsize=10)
            else:
                t, y = ct.step_response(self.C / (1 + sys_ol), T=t_vec)
                self.ax_main.plot(t, y, linewidth=2, label="Active (Current)")
                self.ax_main.set_title("Control Effort", fontsize=10)
            
            # Styles and line options to distinguish saved designs
            styles = ['--', '-.', ':', ':']
            for idx, c_saved in enumerate(self.saved_designs):
                if self.flag_disp_C[idx] == 1:
                    ol_saved, cl_saved = self.get_active_system_eval(c_saved)
                    style = styles[idx % len(styles)]
                    if choice == 0:
                        ts, ys = ct.step_response(cl_saved, T=t_vec)
                        self.ax_main.plot(ts, ys, style, linewidth=1.5, label=f"Design {idx+1}")
                    elif choice == 1:
                        ts, ys = ct.impulse_response(cl_saved, T=t_vec)
                        self.ax_main.plot(ts, ys, style, linewidth=1.5, label=f"Design {idx+1}")
                    else:
                        ts, ys = ct.step_response(c_saved / (1 + ol_saved), T=t_vec)
                        self.ax_main.plot(ts, ys, style, linewidth=1.5, label=f"Design {idx+1}")
                
            self.ax_main.legend(loc="best", fontsize=8, frameon=True)
            self.ax_main.grid(True, linestyle='--', alpha=0.7)
            self.ax_main.set_xlabel("Time (seconds)", fontsize=9)
            self.ax_main.set_ylabel("Amplitude", fontsize=9)
            self.fig_main.tight_layout()
            self.canvas_main.draw()
        except Exception as e:
            print(f"Main plot error: {e}")

        try:
            self.ax_rl.clear()
            ct.root_locus(sys_ol, ax=self.ax_rl)
            self.ax_rl.set_title(f"Root Locus Diagram (Config {self.architecture_idx})", fontsize=10)
            
            for z_val in [0.2, 0.4, 0.6, 0.8]:
                r_vals = np.linspace(0, 50, 100)
                omega_n = r_vals
                real_part = -z_val * omega_n
                imag_part = omega_n * np.sqrt(1 - z_val**2)
                self.ax_rl.plot(real_part, imag_part, ':', color='gray', alpha=0.5)
                self.ax_rl.plot(real_part, -imag_part, ':', color='gray', alpha=0.5)

            for w_val in [2, 5, 10, 20]:
                circ = plt.Circle((0, 0), w_val, color='gray', fill=False, linestyle=':', alpha=0.4)
                self.ax_rl.add_patch(circ)

            self.ax_rl.grid(True, linestyle='--', alpha=0.7)
            self.canvas_rl.draw()
        except Exception as e:
            print(f"Root locus error: {e}")

        try:
            self.plot_bode_complete(sys_ol, self.ax_bode_m, self.ax_bode_p, f"Open-Loop Bode Diagram (Config {self.architecture_idx})", is_grid=False)
            self.fig_bode.tight_layout()
            self.canvas_bode.draw()
        except Exception as e:
            print(f"Bode error: {e}")

        try:
            self.draw_octave_style_nyquist(sys_ol, self.ax_nyq, f"Nyquist Diagram (Config {self.architecture_idx})")
            self.fig_nyq.tight_layout()
            self.canvas_nyq.draw()
        except Exception as e:
            print(f"Nyquist error: {e}")

        try:
            self.ax_all_rl.clear()
            self.ax_all_ny.clear()
            
            ct.root_locus(sys_ol, ax=self.ax_all_rl)
            self.ax_all_rl.set_title("Root Locus", fontsize=9)
            self.ax_all_rl.grid(True, linestyle='--', alpha=0.7)
            
            self.draw_octave_style_nyquist(sys_ol, self.ax_all_ny, "Nyquist")
            self.plot_bode_complete(sys_ol, self.ax_all_bm, self.ax_all_bp, "Bode Mag / Phase", is_grid=True)
            
            self.fig_all.tight_layout()
            self.canvas_all.draw()
        except Exception as e:
            print(f"All plots grid error: {e}")
            
        self.lbl_status.config(text=f"Architecture: Config {self.architecture_idx} | Status: Ready")
        
    def make_active_design(self):
        sel = self.listbox_designs.curselection()
        if not sel:
            messagebox.showwarning("Warning", "Please select a saved design from the list to make active.")
            return
        idx = sel[0]
        sys_saved = self.saved_designs[idx]
        
        try:
            # Extract poles and zeros using python-control functions
            p = ct.poles(sys_saved)
            z = ct.zeros(sys_saved)
            
            self.C_zeros = [self.clean_val(val) for val in np.atleast_1d(z)]
            self.C_poles = [self.clean_val(val) for val in np.atleast_1d(p)]
            
            # Extract gain from transfer function data
            num, den = ct.tfdata(sys_saved)
            num_arr = np.array(num[0][0], dtype=float)
            den_arr = np.array(den[0][0], dtype=float)
            
            num_lead = num_arr[0] if len(num_arr) > 0 else 1.0
            den_lead = den_arr[0] if len(den_arr) > 0 else 1.0
            self.C_gain = float(num_lead / den_lead) if den_lead != 0 else 1.0
                
            # Update active compensator C
            self.C = ct.zpk(self.C_zeros, self.C_poles, self.C_gain)
            self.active_comp_name = "C"
            self.edit_sys = self.C
            
            # Update slider and gain label
            if self.C_gain > self.slider_gain.cget("to"):
                self.slider_gain.configure(to=max(2000, self.C_gain * 2))
            self.slider_gain.set(self.C_gain)
            self.lbl_gain_val.config(text=f"{self.C_gain:.2f}")
            
            # Refresh all plots and the Edit Controller table/equations
            self.refresh_all()
            messagebox.showinfo("Success", f"{self.saved_designs_meta[idx]} is now the active controller!")
        except Exception as e:
            messagebox.showerror("Error", f"Could not load design into active compensator:\n{e}")            

class ToplevelArchitecture(tk.Toplevel):
    def __init__(self, parent_app):
        super().__init__(parent_app.root)
        self.parent_app = parent_app
        self.title("Edit Architecture")
        self.geometry("980x740")
        self.grab_set()

        # Left Frame: Configurations Radio Buttons
        left_frame = ttk.LabelFrame(self, text="Configurations")
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)

        self.arch_var = tk.IntVar(value=self.parent_app.architecture_idx)

        for i in range(1, 7):
            ttk.Radiobutton(left_frame, text=f"Config {i}", variable=self.arch_var, value=i, command=self.on_arch_changed).pack(anchor="w", padx=10, pady=12)

        # Right Frame: Diagram preview + Blocks Table
        right_frame = ttk.Frame(self)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Diagram Preview Canvas Frame
        diag_frame = ttk.LabelFrame(right_frame, text="Architecture Diagram Preview")
        diag_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.fig_arch, self.ax_arch = plt.subplots(figsize=(6, 2.6))
        self.fig_arch.patch.set_facecolor('white')
        self.canvas_arch = FigureCanvasTkAgg(self.fig_arch, master=diag_frame)
        self.canvas_arch.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Blocks Table Frame
        blocks_frame = ttk.LabelFrame(right_frame, text="Blocks")
        blocks_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=5)

        columns = ("id", "name", "value")
        self.tree_blocks = ttk.Treeview(blocks_frame, columns=columns, show="headings", height=6)
        self.tree_blocks.heading("id", text="Identifier")
        self.tree_blocks.heading("name", text="Block Name")
        self.tree_blocks.heading("value", text="Value")
        self.tree_blocks.column("id", width=90, anchor="w")
        self.tree_blocks.column("name", width=180, anchor="w")
        self.tree_blocks.column("value", width=360, anchor="w")
        self.tree_blocks.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        edit_frame = ttk.Frame(right_frame)
        edit_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        ttk.Label(edit_frame, text="Edit Value:").pack(side=tk.LEFT, padx=5)
        self.ent_block_val = ttk.Entry(edit_frame, width=45)
        self.ent_block_val.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        self.block_map = {}
        self.tree_blocks.bind("<<TreeviewSelect>>", self.on_tree_select)

        btn_frame = ttk.Frame(right_frame)
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=5, pady=10)

        ttk.Button(btn_frame, text="OK", command=self.on_ok).pack(side=tk.RIGHT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=self.destroy).pack(side=tk.RIGHT, padx=5)

        self.update_block_fields()
        self.load_and_display_image(self.arch_var.get())

    def on_arch_changed(self):
        cfg = self.arch_var.get()
        self.update_block_fields()
        self.load_and_display_image(cfg)

    def load_and_display_image(self, cfg):
        self.ax_arch.clear()
        self.ax_arch.axis('off')
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        img_path = os.path.join(script_dir, "img", f"Config{cfg}.png")
        
        if os.path.exists(img_path):
            try:
                img_data = plt.imread(img_path)
                self.ax_arch.imshow(img_data)
                self.ax_arch.set_title(f"Config {cfg} Architecture", fontsize=10, fontweight='bold')
            except Exception as e:
                self.ax_arch.text(0.5, 0.5, f"Error loading image:\n{e}", ha='center', va='center', fontsize=9)
        else:
            self.ax_arch.text(0.5, 0.5, f"Image not found:\n{img_path}\n(Please ensure 'img/Config{cfg}.png' exists)", 
                              ha='center', va='center', fontsize=9, wrap=True)
            
        self.fig_arch.tight_layout()
        self.canvas_arch.draw()

    def update_block_fields(self):
        self.tree_blocks.delete(*self.tree_blocks.get_children())
        self.block_map.clear()
        cfg = self.arch_var.get()

        if cfg in [1, 2, 3]:
            items = [
                ("F", "Prefilter F(s)", self.parent_app.get_sys_string("F")),
                ("C", "Controller C(s)", self.parent_app.get_sys_string("C")),
                ("G", "Plant G(s)", self.parent_app.get_sys_string("G")),
                ("H", "Sensor H(s)", self.parent_app.get_sys_string("H"))
            ]
        elif cfg == 4:
            items = [
                ("C1", "Controller C1(s)", self.parent_app.get_sys_string("C1")),
                ("C2", "Controller C2(s)", self.parent_app.get_sys_string("C2")),
                ("G", "Plant G(s)", self.parent_app.get_sys_string("G")),
                ("H", "Sensor H(s)", self.parent_app.get_sys_string("H"))
            ]
        elif cfg == 5:
            items = [
                ("F", "Prefilter F(s)", self.parent_app.get_sys_string("F")),
                ("C", "Controller C(s)", self.parent_app.get_sys_string("C")),
                ("G1", "Plant G1(s)", self.parent_app.get_sys_string("G1")),
                ("G2", "Plant G2(s)", self.parent_app.get_sys_string("G2")),
                ("Gd", "Disturbance Gd(s)", self.parent_app.get_sys_string("Gd"))
            ]
        elif cfg == 6:
            items = [
                ("F", "Prefilter F(s)", self.parent_app.get_sys_string("F")),
                ("C1", "Controller C1(s)", self.parent_app.get_sys_string("C1")),
                ("C2", "Controller C2(s)", self.parent_app.get_sys_string("C2")),
                ("G1", "Plant G1(s)", self.parent_app.get_sys_string("G1")),
                ("G2", "Plant G2(s)", self.parent_app.get_sys_string("G2")),
                ("H1", "Sensor H1(s)", self.parent_app.get_sys_string("H1")),
                ("H2", "Sensor H2(s)", self.parent_app.get_sys_string("H2"))
            ]
        else:
            items = []

        for item in items:
            item_id = self.tree_blocks.insert("", tk.END, values=item)
            self.block_map[item_id] = item[0]

        while len(items) < 7:
            self.tree_blocks.insert("", tk.END, values=("--", "--", "--"))
            items.append(None)

    def on_tree_select(self, event):
        sel = self.tree_blocks.selection()
        if sel:
            vals = self.tree_blocks.item(sel[0], "values")
            if vals and vals[0] != "--":
                self.ent_block_val.delete(0, tk.END)
                self.ent_block_val.insert(0, vals[2])

    def on_ok(self):
        cfg = self.arch_var.get()
        self.parent_app.architecture_idx = cfg
        s = ct.tf('s')

        for item_id in self.tree_blocks.get_children():
            vals = self.tree_blocks.item(item_id, "values")
            bid = vals[0]
            if bid != "--":
                sel = self.tree_blocks.selection()
                if sel and sel[0] == item_id:
                    new_val = self.ent_block_val.get()
                else:
                    new_val = vals[2]

                try:
                    parsed = eval(new_val)
                    if isinstance(parsed, (int, float, complex)):
                        parsed = ct.tf([parsed], [1])
                    setattr(self.parent_app, bid, parsed)
                    setattr(self.parent_app, f"{bid}_str", new_val)
                    if bid == "G" and hasattr(self.parent_app, 'ent_plant'):
                        self.parent_app.ent_plant.delete(0, tk.END)
                        self.parent_app.ent_plant.insert(0, new_val)
                except Exception as e:
                    messagebox.showerror("Error", f"Invalid syntax for block {bid}: {e}", parent=self)
                    return

        self.destroy()
        self.parent_app.refresh_all()   

if __name__ == "__main__":
    root = tk.Tk()
    app = SISOToolApp(root)
    root.mainloop()
