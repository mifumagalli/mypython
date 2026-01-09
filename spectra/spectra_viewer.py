import sys
import numpy as np
import matplotlib.pyplot as plt
import os

def load_spectrum(filename):
    """Loads spectrum from an ASCII file (wavelength, flux, error)."""
    try:
        data = np.loadtxt(filename)
        if data.shape[1] < 2:
            raise ValueError("File must have at least 2 columns (wavelength, flux)")
        wavelength = data[:, 0]
        flux = data[:, 1]
        error = data[:, 2] if data.shape[1] > 2 else np.zeros_like(flux)
        return wavelength, flux, error
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

from matplotlib.widgets import Button, TextBox, CheckButtons, RadioButtons
from linetools.lists.linelist import LineList
from astropy.stats import sigma_clipped_stats
import warnings

# Suppress linetools warnings if any
warnings.filterwarnings("ignore")

class SpectraViewer:
    def __init__(self, filename):
        self.filename = filename
        self.wavelength, self.flux, self.error = load_spectrum(filename)
        self.linelist = LineList('Strong')
        self.lines_visible = False
        self.current_z = 0.0
        self.vlines = []

        # EW State: Key = (name, wrest.value), Value = {'selected': bool, 'custom_window': (min, max) or None}
        self.line_ew_state = {}
        # Setup Figure and Gridspec
        self.fig = plt.figure(figsize=(12, 6))
        # Main plot takes 80%, Control panel takes 20%
        self.gs = self.fig.add_gridspec(1, 2, width_ratios=[4, 1])
        
        self.ax = self.fig.add_subplot(self.gs[0])
        self.line, = self.ax.plot(self.wavelength, self.flux, 'b-', label='Flux')
        self.error_line, = self.ax.plot(self.wavelength, self.error, 'orange', label='Error', alpha=0.7)
        
        self.ax.set_xlabel('Wavelength')
        self.ax.set_ylabel('Flux')
        self.ax.set_title(f'Spectrum: {filename}')
        self.ax.legend()
        
        # Helper for initial view
        self.ax.autoscale(enable=True, axis='both', tight=True)
        self.ax.relim()
        self.ax.autoscale_view()

        # Control Panel
        self.panel_ax = self.fig.add_subplot(self.gs[1])
        self.panel_ax.axis('off')

        # --- Spectrum Controls ---
        # Use figure coords to align with widgets (Left ~0.82)
        self.fig.text(0.82, 0.94, "Spectrum Controls", fontweight='bold', fontsize=10)
        line1 = plt.Line2D([0.82, 0.97], [0.92, 0.92], transform=self.fig.transFigure, color='k', linewidth=1)
        self.fig.add_artist(line1)
        
        # Checkbox for showing lines
        # y: 0.85 - 0.90
        ax_check = self.fig.add_axes([0.82, 0.85, 0.15, 0.05])
        self.check = CheckButtons(ax_check, ['Show Lines'], [self.lines_visible])
        self.check.on_clicked(self.toggle_lines)
        
        # Textbox for Redshift
        # y: 0.75 - 0.80
        ax_text = self.fig.add_axes([0.82, 0.75, 0.15, 0.05])
        self.text_box = TextBox(ax_text, 'Redshift: ', initial=str(self.current_z))
        self.text_box.on_submit(self.submit_redshift)

        # --- Continuum Fit ---
        # Header at 0.65 -> Move up slightly to 0.68? Or keep 0.65 but move widgets down.
        # Let's clean up spacing.
        self.fig.text(0.82, 0.66, "Continuum Fit", fontweight='bold', fontsize=10)
        line2 = plt.Line2D([0.82, 0.97], [0.64, 0.64], transform=self.fig.transFigure, color='k', linewidth=1)
        self.fig.add_artist(line2)

        # Continuum Mode
        self.continuum_mode = False
        self.continuum_anchors = [] 
        self.continuum_plots = [] 
        self.continuum_fit_line = None
        self.poly_order = 1
        
        # Frozen Segments
        self.frozen_segments = [] 
        self.freeze_start = None
        self.freeze_guide = None 
        self.model_accepted = False

        # Checkbox for Continuum Mode
        # y: 0.55 -> Move to 0.57? 
        # Header at 0.66. Line at 0.64. Widget should be below 0.64.
        # Widget y=0.55, h=0.05 -> Top 0.60. Spacing 0.04. Good.
        ax_check_cont = self.fig.add_axes([0.82, 0.55, 0.15, 0.05])
        self.check_cont = CheckButtons(ax_check_cont, ['Continuum'], [self.continuum_mode])
        self.check_cont.on_clicked(self.toggle_continuum_mode)
        
        # Textbox for Poly Order
        # y: 0.45 -> Move to 0.47? Standard spacing of 0.08?
        # 0.55 - 0.08 = 0.47.
        # Let's keep 0.45 (Gap 0.10).
        ax_text_poly = self.fig.add_axes([0.82, 0.45, 0.15, 0.05])
        self.text_poly = TextBox(ax_text_poly, 'Poly Order: ', initial=str(self.poly_order))
        self.text_poly.on_submit(self.submit_poly_order)

        # Auto Continuum Button
        # y: 0.25? No, putting it between Poly Order (0.45) and Save (0.35)?
        # Or below Save?
        # User said "in continuum part".
        # Let's put it above Save.
        # Save is at 0.35. Let's move Save down or put Auto at 0.25?
        # Actually ample space below 0.35.
        # Let's putting Auto at 0.25.
        
        # Auto Continuum Button
        # Moved up from 0.25 to 0.37 to close gap from PolyOrder (0.45)
        ax_auto = self.fig.add_axes([0.82, 0.37, 0.15, 0.05])
        self.btn_auto = Button(ax_auto, 'Auto Cont')
        self.btn_auto.on_clicked(self.auto_continuum)
        
        # Save Fit Button
        # Moved up from 0.18 to 0.30
        ax_save = self.fig.add_axes([0.82, 0.30, 0.15, 0.05])
        self.btn_save = Button(ax_save, 'Accept fit to save')
        self.btn_save.on_clicked(self.save_fit)
        self.btn_save.label.set_fontsize(8)
        self.btn_save.color = '0.9' 
        self.btn_save.hovercolor = '0.9'

        # --- Equivalent Width ---
        # Moved up from 0.15 to 0.25
        self.fig.text(0.82, 0.25, "Equivalent Width", fontweight='bold', fontsize=10)
        line3 = plt.Line2D([0.82, 0.97], [0.23, 0.23], transform=self.fig.transFigure, color='k', linewidth=1)
        self.fig.add_artist(line3)
        
        # Checkbox for Linear Fallback
        # Moved up from 0.10 to 0.19
        ax_check_ew = self.fig.add_axes([0.82, 0.19, 0.15, 0.03]) # Smaller height
        self.linear_cont_fallback = False # Default Off
        self.check_ew = CheckButtons(ax_check_ew, ['Linear Cont'], [self.linear_cont_fallback])
        self.check_ew.on_clicked(self.toggle_linear_fallback)
        
        # Ew Labels (using figure text for simple display)
        # Moved up from 0.07/0.05 to 0.15/0.13
        self.text_ew_obs = self.fig.text(0.82, 0.15, "Obs: --", fontsize=9)
        self.text_ew_rest = self.fig.text(0.82, 0.13, "Rest: --", fontsize=9)
        
        self.ew_start = None
        self.ew_guide = None
        self.ew_poly = None # Store plotted area polygon
        self.ew_obs_value = None # Store last calculated Obs EW
        
        # Help Button
        ax_help = self.fig.add_axes([0.82, 0.02, 0.15, 0.05])
        self.btn_help = Button(ax_help, 'Help')
        self.btn_help.on_clicked(self.show_help)

        self.connect()

        print("Controls:")
        print("  i / o : Zoom in / out")
        print("  [ / ] : Pan left / right")
        print("  b     : Set bottom Y limit to cursor position")
        print("  t     : Set top Y limit to cursor position")
        print("  w     : Reset window to default view")
        print("  m     : Identify line at cursor")
        print("  z     : Open Line Check Window (3x2 grid of lines)")
        print("  e     : Equivalent Width (press twice for range)")
        print("  l/L   : Toggle log scales on y/x axis")
        print("  q     : Quit GUI")
        print("\nContinuum Mode Controls (when enabled):")
        print("  x     : Add anchor point")
        print("  d     : Delete closest anchor point")
        print("  c     : Compute/Update fit")
        print("  Space : Start/End freeze range (press twice)")
        print("  u     : Unfreeze segment at cursor")
        print("  a     : Accept All (freeze remaining fit)")
        
        plt.show()

    def connect(self):
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.ax.callbacks.connect('xlim_changed', self.on_xlim_changed)

    def on_xlim_changed(self, event_ax):
        if self.lines_visible:
            self.draw_lines()

    def toggle_lines(self, label):
        self.lines_visible = not self.lines_visible
        self.draw_lines()

    def toggle_continuum_mode(self, label):
        self.continuum_mode = not self.continuum_mode
        if not self.continuum_mode:
            # Clear anchors and fit if exiting mode
            # Or keep them visual? Usually exiting mode clears temporary visuals or keeps them?
            # User might want to toggle mode to interact with other things and come back.
            # But let's say visuals persist until cleared explicitly?
            # Requirement doesn't specify. I'll persist data but maybe hide if desired?
            # For now, just toggle state.
            pass

    def toggle_linear_fallback(self, label):
        self.linear_cont_fallback = not self.linear_cont_fallback
        print(f"Linear Continuum Fallback: {self.linear_cont_fallback}")

    def submit_poly_order(self, text):
        try:
            self.poly_order = int(text)
        except ValueError:
            print("Invalid order value")

    def submit_redshift(self, text):
        try:
            self.current_z = float(text)
            self.draw_lines()
            
            # Update Rest EW if available
            if self.ew_obs_value is not None:
                ew_rest = self.ew_obs_value / (1 + self.current_z)
                self.text_ew_rest.set_text(f"Rest: {ew_rest:.3f} A")
                self.fig.canvas.draw_idle()
                
        except ValueError:
            print("Invalid redshift value")

    def draw_lines(self):
        # Remove existing lines
        for line in self.vlines:
            line.remove()
        self.vlines = []
        
        if not self.lines_visible:
            self.fig.canvas.draw_idle()
            return

        # Get wavelength limits to optimize drawing
        xlim = self.ax.get_xlim()
        
        # Calculate observed wavelengths
        # z = (obs / rest) - 1  => obs = rest * (1 + z)
        observed_waves = self.linelist.wrest.value * (1 + self.current_z)
        
        # Filter lines within view (plus some buffer)
        mask = (observed_waves > xlim[0]) & (observed_waves < xlim[1])
        visible_waves = observed_waves[mask]
        visible_names = self.linelist.name[mask]

        # Draw lines
        ymin, ymax = self.ax.get_ylim()
        for wv, name in zip(visible_waves, visible_names):
            vline = self.ax.axvline(wv, color='k', linestyle='--', alpha=0.5)
            self.vlines.append(vline)
            
            # Add label near top
            # We use transform=self.ax.get_xaxis_transform() so y is in axes coordinates (0-1)
            # But x is data coordinates.
            text = self.ax.text(wv, 0.95, name, rotation=90, verticalalignment='top', 
                                horizontalalignment='right', transform=self.ax.get_xaxis_transform(),
                                fontsize=8, clip_on=True)
            self.vlines.append(text)
            
        self.fig.canvas.draw_idle()

    def identify_line(self, cursor_x):
        # We need to ask the user which line this corresponds to.
        # Sort lines by rest wavelength
        sorted_indices = self.linelist.wrest.argsort()
        sorted_names = self.linelist.name[sorted_indices]
        sorted_wrest = self.linelist.wrest[sorted_indices]
        
        # Just show all strong lines for selection
        labels = [f"{name} ({wrest.value:.1f})" for name, wrest in zip(sorted_names, sorted_wrest)]

        # Calculate chunks for columns
        N_per_col = 15
        COLS_PER_PAGE = 3
        
        # Adjust figure size for fixed number of columns
        width = 2.5 * COLS_PER_PAGE
        height = min(6, 0.4 * N_per_col) + 1.0 # Add space for buttons
        
        id_fig = plt.figure(figsize=(width, height))
        id_fig.canvas.manager.set_window_title("Select Line")
        
        self.current_page = 0
        self.radio_buttons = [] 
        self.page_buttons = [] # Store prev/next buttons
        
        def draw_page():
            # Clear axes
            id_fig.clf()
            self.radio_buttons = []
            self.page_buttons = []
            
            # Recalculate layout
            start_x = 0.05
            col_width = 0.9 / COLS_PER_PAGE
            
            # Determine range for this page
            page_start_idx = self.current_page * N_per_col * COLS_PER_PAGE
            
            # Draw columns
            for i in range(COLS_PER_PAGE):
                col_offset = i * N_per_col
                start = page_start_idx + col_offset
                if start >= len(labels): break
                
                end = min(start + N_per_col, len(labels))
                chunk_labels = labels[start:end]
                
                # Position: [left, bottom, width, height]
                # Reserve bottom 0.1 for buttons
                ax_col = id_fig.add_axes([start_x + i*col_width, 0.15, col_width, 0.8])
                ax_col.axis('off')

                radio = RadioButtons(ax_col, chunk_labels, active=0)
                self.radio_buttons.append(radio)
                
                # Closure to capture correct labels
                def make_callback(offset):
                    def on_select(label):
                        if label not in labels: return 
                        idx = labels.index(label)
                        
                        rest_wave = sorted_wrest[idx].value
                        new_z = (cursor_x / rest_wave) - 1
                        
                        print(f"Selected {sorted_names[idx]} at rest {rest_wave:.2f}. Obs: {cursor_x:.2f}. New z: {new_z:.4f}")
                        
                        self.current_z = new_z
                        self.text_box.set_val(f"{new_z:.4f}")
                        
                        if not self.lines_visible:
                            self.check.set_active(0)
                        
                        plt.close(id_fig)
                    return on_select

                radio.on_clicked(make_callback(start))

            # Add navigation buttons
            ax_prev = id_fig.add_axes([0.1, 0.02, 0.3, 0.08])
            b_prev = Button(ax_prev, 'Previous')
            
            ax_next = id_fig.add_axes([0.6, 0.02, 0.3, 0.08])
            b_next = Button(ax_next, 'Next')
            
            self.page_buttons.extend([b_prev, b_next])
            
            def prev_clicked(event):
                if self.current_page > 0:
                    self.current_page -= 1
                    draw_page()
            
            def next_clicked(event):
                # Check if next page exists
                if (self.current_page + 1) * N_per_col * COLS_PER_PAGE < len(labels):
                    self.current_page += 1
                    draw_page()

            b_prev.on_clicked(prev_clicked)
            b_next.on_clicked(next_clicked)
            
            id_fig.canvas.draw()

        draw_page()
        plt.show()

    def show_line_check_window(self):
        # Determine observed range
        min_obs = np.min(self.wavelength)
        max_obs = np.max(self.wavelength)
        
        # Calculate observed wavelengths for all lines
        # linelist.wrest is a Quantity array (Astropy)
        all_obs_waves = self.linelist.wrest.value * (1 + self.current_z)
        
        # Filter
        in_range_mask = (all_obs_waves >= min_obs) & (all_obs_waves <= max_obs)
        
        if not np.any(in_range_mask):
             print(f"No lines found in range {min_obs:.1f} - {max_obs:.1f} A for z={self.current_z:.4f}")
             return

        # Get filtered subset
        filtered_names = self.linelist.name[in_range_mask]
        filtered_wrest = self.linelist.wrest[in_range_mask]
        
        # Sort lines by rest wavelength
        sorted_indices = filtered_wrest.argsort()
        sorted_names = filtered_names[sorted_indices]
        sorted_wrest = filtered_wrest[sorted_indices]
        
        # 3x2 grid
        COLS = 3
        ROWS = 2
        N_PER_PAGE = COLS * ROWS
        
        # Figure setup
        check_fig = plt.figure(figsize=(15, 8))
        check_fig.canvas.manager.set_window_title(f"Line Check (z={self.current_z:.4f})")
        
        self.check_current_page = 0
        self.check_buttons = []
        
        def draw_check_page():
            # Clear axes and map
            check_fig.clf()
            self.check_buttons = []
            self.check_axes_map = {}
            
            # Subplots
            # We want to use add_subplot usually, but let's manage them carefully to allow buttons at bottom
            # Reserve bottom 0.1 for buttons
            gs = check_fig.add_gridspec(ROWS, COLS, bottom=0.15, top=0.95, hspace=0.4, wspace=0.3)
            
            start_idx = self.check_current_page * N_PER_PAGE
            
            for i in range(N_PER_PAGE):
                idx = start_idx + i
                if idx >= len(sorted_names): break
                
                rest_wave = sorted_wrest[idx].value
                name = sorted_names[idx]
                
                obs_wave = rest_wave * (1 + self.current_z)
                
                # +/- 1000 km/s window
                # dv/c = dlambda/lambda => dlambda = lambda * (dv/c)
                c_kms = 299792.458
                delta_wave = obs_wave * (1000.0 / c_kms)
                
                w_min = obs_wave - delta_wave
                w_max = obs_wave + delta_wave
                
                # Get data in range
                mask = (self.wavelength >= w_min) & (self.wavelength <= w_max)
                
                row = i // COLS
                col = i % COLS
                ax = check_fig.add_subplot(gs[row, col])
                
                # Store map
                self.check_axes_map[ax] = (name, rest_wave)
                
                if np.any(mask):
                    ax.plot(self.wavelength[mask], self.flux[mask], 'b-', linewidth=1)
                    ax.plot(self.wavelength[mask], self.error[mask], 'orange', alpha=0.5, linewidth=1)
                
                # Center line
                ax.axvline(obs_wave, color='k', linestyle='--')
                
                ax.set_title(f"{name}\nRest: {rest_wave:.1f}, Obs: {obs_wave:.1f}", fontsize=10)
                ax.set_xlim(w_min, w_max)
                # Auto ylim if data exists
                if np.any(mask):
                    y_seg = self.flux[mask]
                    # Robust scaling?
                    # simple min max for now
                    if len(y_seg) > 0:
                         ymin, ymax = np.min(y_seg), np.max(y_seg)
                         margin = (ymax - ymin) * 0.1 if ymax != ymin else 1.0
                         ax.set_ylim(ymin - margin, ymax + margin)

            # Navigation Buttons
            # Prev
            if self.check_current_page > 0:
                ax_prev = check_fig.add_axes([0.1, 0.05, 0.1, 0.05])
                b_prev = Button(ax_prev, 'Previous')
                b_prev.on_clicked(lambda e: change_page(-1))
                self.check_buttons.append(b_prev)
            
            # Next
            if (self.check_current_page + 1) * N_PER_PAGE < len(sorted_names):
                ax_next = check_fig.add_axes([0.8, 0.05, 0.1, 0.05])
                b_next = Button(ax_next, 'Next')
                b_next.on_clicked(lambda e: change_page(1))
                self.check_buttons.append(b_next)
                
            # Measure EW Button
            ax_measure = check_fig.add_axes([0.65, 0.05, 0.12, 0.05])
            b_measure = Button(ax_measure, 'Measure EW')
            b_measure.on_clicked(self.measure_ew_batch)
            self.check_buttons.append(b_measure)
                
            # Page Info
            check_fig.text(0.5, 0.075, f"Page {self.check_current_page + 1} / {int(np.ceil(len(sorted_names)/N_PER_PAGE))}", 
                           horizontalalignment='center', fontsize=12)

            # Usage Hint
            hint_text = "Controls: 'u' Select | 'd' Deselect | 'e'+'e' Custom Window"
            check_fig.text(0.5, 0.02, hint_text, horizontalalignment='center', fontsize=10, color='blue')
            
            check_fig.canvas.draw()
            
        def change_page(delta):
            self.check_current_page += delta
            draw_check_page()
            
        # Connect key press for this window
        check_fig.canvas.mpl_connect('key_press_event', self.on_check_key)
        
        # Store for key handler access
        self.check_fig = check_fig
        self.sorted_names_cache = sorted_names
        self.sorted_wrest_cache = sorted_wrest
        self.check_axes_cache = {} # Map (row, col) to data index or similar? 
        # Actually better to map axes to line info.
        
        draw_check_page()
        plt.show()

    def on_check_key(self, event):
        if not event.inaxes: return
        
        # Identify which sub-plot (line) we are in
        # We need a way to map axes to the specific line index.
        # Let's simple iterate over axes in the figure?
        # Or store metadata on axes?
        # Matplotlib axes don't easily store custom data, but we can check equality.
        
        # Re-derive index from logical layout is tricky if we don't store it.
        # Let's modify draw_check_page to store a map: self.check_axes_map[ax] = (name, wrest)
        
        if not hasattr(self, 'check_axes_map') or event.inaxes not in self.check_axes_map:
             return
             
        name, wrest_val = self.check_axes_map[event.inaxes]
        key_id = (name, wrest_val)
        
        # Ensure state entry exists
        if key_id not in self.line_ew_state:
             self.line_ew_state[key_id] = {'selected': False, 'custom_window': None}
        
        state = self.line_ew_state[key_id]
        
        if event.key == 'u':
             state['selected'] = True
             self.update_check_panel_title(event.inaxes, name, wrest_val, state)
             
        elif event.key == 'd':
             state['selected'] = False
             state['custom_window'] = None # Reset custom window on deselect? Or keep?
             # User says "excluded from list". 
             self.update_check_panel_title(event.inaxes, name, wrest_val, state)
             
        elif event.key == 'e':
             # Custom window logic
             # We need a temporary state: "waiting for second 'e'"
             if 'ew_wait_start' not in state:
                 state['ew_wait_start'] = event.xdata
                 # Draw guide
                 guide = event.inaxes.axvline(event.xdata, color='green', linestyle='--')
                 state['guide_artist'] = guide
                 event.inaxes.set_title(f"Define End...", color='blue') # Feedback
                 event.canvas.draw_idle()
             else:
                 # End
                 start = min(state['ew_wait_start'], event.xdata)
                 end = max(state['ew_wait_start'], event.xdata)
                 state['custom_window'] = (start, end)
                 
                 # Cleanup
                 if 'guide_artist' in state:
                     state['guide_artist'].remove()
                     del state['guide_artist']
                 del state['ew_wait_start']
                 
                 # Auto-select if defining window? Probably yes.
                 state['selected'] = True
                 
                 self.draw_check_panel_extras(event.inaxes, state)
                 self.update_check_panel_title(event.inaxes, name, wrest_val, state)

    def update_check_panel_title(self, ax, name, wrest, state):
        obs = wrest * (1 + self.current_z)
        title_text = f"{name}\nRest: {wrest:.1f}, Obs: {obs:.1f}"
        
        color = 'black'
        if state.get('selected'):
             color = 'green'
        
        # Explicit check for FALSE vs not present? 
        # Logic: 'u' -> True (included). 'd' -> False (excluded).
        # Default? If unknown, black.
        
        if state.get('selected') is False: # Explicitly False
             color = 'red'
             
        ax.set_title(title_text, color=color, fontsize=10)
        self.check_fig.canvas.draw_idle()

    def draw_check_panel_extras(self, ax, state):
        # Clear previous extras (shaded regions)
        # We need to track them. Let's clear specific artists if we stored them, or just rely on redraw?
        # Redrawing the whole page is expensive.
        # Let's try to manage just the patch.
        
        if 'window_patch' in state and state['window_patch']:
             try:
                 state['window_patch'].remove()
             except: pass
             state['window_patch'] = None
             
        if state.get('custom_window'):
             w_min, w_max = state['custom_window']
             # Draw shaded region
             patch = ax.axvspan(w_min, w_max, color='green', alpha=0.2)
             state['window_patch'] = patch
             
        # Also maybe default window?
        # Requirement: "Default of +/- 200 km/s".
        # Should we visualize the default if selected but no custom window?
        # User didn't explicitly ask, but it's helpful.
        # "Title turns green". 
        
        self.check_fig.canvas.draw_idle()

    def measure_ew_batch(self, event):
        print("Measuring EW for selected lines...")
        
        results = []
        
        # Iterate over all possible lines or just state?
        # State might only contain lines we interacted with.
        # User might want to 'u' a line.
        # If a line is NOT in state, is it selected? Default NO.
        
        for key_id, state in self.line_ew_state.items():
             if state.get('selected'):
                 name, wrest = key_id
                 
                 # Determine window
                 if state.get('custom_window'):
                     w_start, w_end = state['custom_window']
                     is_custom = True
                 else:
                     # Default +/- 200 km/s
                     obs_wave = wrest * (1 + self.current_z)
                     c_kms = 299792.458
                     delta = obs_wave * (200.0 / c_kms)
                     w_start = obs_wave - delta
                     w_end = obs_wave + delta
                     is_custom = False
                     
                 # Calculation
                 # Reuse logic from handle_ew? Or copy for batch efficiency/customization.
                 
                 # Get data slice
                 mask = (self.wavelength >= w_start) & (self.wavelength <= w_end)
                 if not np.any(mask):
                     print(f"Skipping {name}: No data in window.")
                     continue
                 
                 wave_seg = self.wavelength[mask]
                 flux_seg = self.flux[mask]
                 err_seg = self.error[mask]
                 
                 # Determine Continuum
                 # 1. Global model
                 cont_seg = None
                 if hasattr(self, 'continuum_model') and not np.all(np.isnan(self.continuum_model[mask])):
                     cont_seg = self.continuum_model[mask]
                 elif hasattr(self, 'current_poly') and self.current_poly is not None:
                     cont_seg = self.current_poly(wave_seg)
                 else:
                     # Fallback 1.0
                     cont_seg = np.ones_like(flux_seg)
                     
                 # EW
                 norm_flux = flux_seg / cont_seg
                 integrand = 1.0 - norm_flux
                 
                 ew_obs = np.trapz(integrand, x=wave_seg) # Positive for absorption
                 
                 # Error
                 # Propagate error: sigma_EW^2 = sum( (sigma_flux / Cont * dlambda)^2 ) approximately
                 # assuming independent errors per pixel.
                 # dlambda can be approximated by gradient or diff
                 # Let's use simple mid-point spacing
                 dl = np.gradient(wave_seg)
                 term = (err_seg / cont_seg) * dl
                 ew_err_obs = np.sqrt(np.sum(term**2))
                 
                 ew_rest = ew_obs / (1 + self.current_z)
                 ew_err_rest = ew_err_obs / (1 + self.current_z)
                 
                 obs_center = wrest * (1 + self.current_z)
                 
                 # Columns: Name, RestWave, ObsWave, WinStart, WinEnd, EW_Obs, Err_Obs, EW_Rest, Err_Rest
                 results.append((name, wrest, obs_center, w_start, w_end, ew_obs, ew_err_obs, ew_rest, ew_err_rest))
                 
        if not results:
             print("No lines selected or valid for measurement.")
             return
             
        # Write to file
        fname = os.path.splitext(self.filename)[0] + "_EW.txt"
        
        try:
             with open(fname, 'w') as f:
                 f.write("# Name RestWave ObsWave WinStart WinEnd EW_Obs Err_Obs EW_Rest Err_Rest\n")
                 for res in results:
                     # Format: ensure Name handles spaces? usually name has no spaces or we quote
                     # Assuming name is safe.
                     line_fmt = "{:s} {:.2f} {:.2f} {:.2f} {:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n"
                     f.write(line_fmt.format(*res))
             print(f"Saved measurements for {len(results)} lines to {fname}")
        except Exception as e:
             print(f"Error saving EW file: {e}")

    def on_key(self, event):
        if not event.inaxes and event.key != 'l': # 'l' might be pressed without axes focus? No, needs cursor position.
             if event.key != 'l': return

        # ... existing logic ...

        # Get current limits
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]

        key = event.key
        
        # Zooming
        if key == 'i': # Zoom In
            # Zoom centered on cursor
            scale_factor = 0.5 # Zoom in by 2x
            new_range = x_range * scale_factor
            center = event.xdata
            # Re-center logic: try to keep cursor at same relative position or just center?
            # Requirement: "zoom scale centered on cursor" is usually preferred.
            # Let's simple center on cursor for now as it's intuitive.
            new_min = center - new_range / 2
            new_max = center + new_range / 2
            self.ax.set_xlim(new_min, new_max)
            
        elif key == 'o': # Zoom Out
            scale_factor = 2.0 # Zoom out by 2x
            new_range = x_range * scale_factor
            center = event.xdata
            new_min = center - new_range / 2
            new_max = center + new_range / 2
            self.ax.set_xlim(new_min, new_max)

        # Panning
        elif key == '[': # Pan Left
            pan_amount = x_range * 0.2
            self.ax.set_xlim(xlim[0] - pan_amount, xlim[1] - pan_amount)
            
        elif key == ']': # Pan Right
            pan_amount = x_range * 0.2
            self.ax.set_xlim(xlim[0] + pan_amount, xlim[1] + pan_amount)

        # Vertical Limits
        elif key == 'b': # Set bottom
            self.ax.set_ylim(bottom=event.ydata)
            
        elif key == 't': # Set top
            self.ax.set_ylim(top=event.ydata)

        elif key == 'l':
             self.ax.set_yscale('log' if self.ax.get_yscale() == 'linear' else 'linear')
        elif key == 'L':
             self.ax.set_xscale('log' if self.ax.get_xscale() == 'linear' else 'linear')
        elif key == 'q':
             plt.close('all')
             sys.exit(0)

        # Reset View
        if key == 'w': # Reset to default
            self.ax.autoscale(enable=True, axis='both', tight=True)
            self.ax.relim()
            self.ax.autoscale_view()
            self.draw_lines() # Re-draw lines to match new view if needed? No, lines are absolute coords.
            
        elif key == 'm': # Identify line
             if event.inaxes == self.ax:
                 self.identify_line(event.xdata)

        # Continuum Mode Keys
        if self.continuum_mode and event.inaxes == self.ax:
            if key == 'x':
                self.add_anchor(event.xdata, event.ydata)
            elif key == 'd':
                self.delete_anchor(event.xdata, event.ydata)
            elif key == 'c':
                self.update_fit()
            elif key == ' ':
                self.handle_freeze(event.xdata)
            elif key == 'u':
                self.handle_unfreeze(event.xdata)
            elif key == 'a':
                self.handle_accept_all()
        
        # We should allow 'e' even if not in continuum mode?
        # Requirement says "This mode is activated by hitting 'e'".
        # Implies it might be independent.
        # But control panel section implies it's always there.
        # Let's assume global access or continuum mode? 
        # Usually analysis like EW is general.
        # But my implementation relies on handle_ew being called.
        # Let's move it OUTSIDE continuum_mode block allow it generally or check my previous 'm' (Identify) logic.
        # 'm' is separate.
        # Let's put 'e' at same level as 'm'.
        
        if key == 'e' and event.inaxes == self.ax:
             self.handle_ew(event.xdata)

        if key == 'z':
             self.show_line_check_window()

        self.fig.canvas.draw_idle()

    def handle_freeze(self, x):
        if self.freeze_start is None:
            # Start freezing range
            self.freeze_start = x
            self.freeze_guide = self.ax.axvline(x, color='grey', linestyle='--', alpha=0.5)
            print(f"Freeze start at {x:.2f}. Press Space again to end.")
        else:
            # End freezing range
            start = min(self.freeze_start, x)
            end = max(self.freeze_start, x)
            print(f"Freezing range: {start:.2f} - {end:.2f}")
            
            # Remove guide
            if self.freeze_guide:
                self.freeze_guide.remove()
                self.freeze_guide = None
            self.freeze_start = None
            
            # Capture current fit
            if not self.continuum_fit_line:
                print("No active fit to freeze.")
                return
            
            # The fit line might be on a dense grid (from update_fit before refactor) 
            # OR we should just evaluate the polynomial on the wavelength grid indices corresponding to the range.
            # But the current fit polynomial is needed.
            
            # Let's verify if we have a valid polynomial fit first.
            if not hasattr(self, 'current_poly') or self.current_poly is None:
                 print("No active polynomial fit.")
                 return

            # Identify indices in self.wavelength
            mask_indices = np.where((self.wavelength >= start) & (self.wavelength <= end))[0]
            if len(mask_indices) == 0:
                print("No data points in frozen range.")
                return
            
            # Evaluate fit on these points
            y_values = self.current_poly(self.wavelength[mask_indices])
            
            # Store segment: mask (indices) and y values
            self.frozen_segments.append({
                'mask_idx': mask_indices,
                'y': y_values,
                'range': (start, end)
            })
            
            self.update_continuum_display()

    def handle_unfreeze(self, x):
        # check if x is in any frozen segment range
        to_remove = []
        for i, seg in enumerate(self.frozen_segments):
            if seg['range'][0] <= x <= seg['range'][1]:
                to_remove.append(i)
        
        for i in sorted(to_remove, reverse=True):
            print(f"Unfreezing segment {i}")
            self.frozen_segments.pop(i)
            
        if to_remove:
            self.update_continuum_display()

    def handle_accept_all(self):
        # Freeze remaining non-frozen regions with current fit
        if not hasattr(self, 'current_poly') or self.current_poly is None:
             print("No active fit to accept.")
             return
             
        # Find indices that are NOT covered by any frozen segment
        # We can construct a coverage array
        coverage = np.zeros(len(self.wavelength), dtype=bool)
        for seg in self.frozen_segments:
            coverage[seg['mask_idx']] = True
            
        # Indices to freeze
        gap_indices = np.where(~coverage)[0]
        
        if len(gap_indices) > 0:
            # We might have disjoint gaps, so simple approach: Find contiguous chunks
            diffs = np.diff(gap_indices)
            breaks = np.where(diffs > 1)[0]
            
            starts = np.insert(breaks + 1, 0, 0)
            ends = np.append(breaks, len(gap_indices) - 1)
            
            added_count = 0
            for s, e in zip(starts, ends):
                indices = gap_indices[s:e+1]
                if len(indices) == 0: continue
                
                start_w = self.wavelength[indices[0]]
                end_w = self.wavelength[indices[-1]]
                y_vals = self.current_poly(self.wavelength[indices])
                
                self.frozen_segments.append({
                    'mask_idx': indices,
                    'y': y_vals,
                    'range': (start_w, end_w)
                })
                added_count += 1
                
            print(f"Accepted fit into {added_count} new segments.")
            self.update_continuum_display()

        # Mark as accepted and enable Save
        self.model_accepted = True
        self.btn_save.label.set_text("Save Fit")
        self.btn_save.color = '0.85' # Standard button color
        self.btn_save.hovercolor = '0.95' # Hover effect
        print("Model Accepted. You can now save the fit.")
        self.fig.canvas.draw_idle()

    def update_continuum_display(self):
        # Combine all frozen segments
        # We need to compute Mean Y at each wavelength point
        
        if not self.frozen_segments:
             # Just clear any existing combined line
             if hasattr(self, 'frozen_line_obj') and self.frozen_line_obj:
                 self.frozen_line_obj.remove()
                 self.frozen_line_obj = None
             self.fig.canvas.draw_idle()
             return

        sum_y = np.zeros(len(self.wavelength))
        count_y = np.zeros(len(self.wavelength))
        
        for seg in self.frozen_segments:
            idx = seg['mask_idx']
            sum_y[idx] += seg['y']
            count_y[idx] += 1
            
        # Compute mean where count > 0
        mask_valid = count_y > 0
        self.continuum_model = np.zeros_like(sum_y) # Store for saving
        self.continuum_model[:] = np.nan
        
        self.continuum_model[mask_valid] = sum_y[mask_valid] / count_y[mask_valid]
        
        # Plot
        # Remove old line
        if hasattr(self, 'frozen_line_obj') and self.frozen_line_obj:
             self.frozen_line_obj.remove()
        
        self.frozen_line_obj, = self.ax.plot(self.wavelength, self.continuum_model, color='grey', linewidth=2, linestyle='-')
        self.fig.canvas.draw_idle()
        
    def save_fit(self, event):
        if not self.model_accepted:
            print("Please press 'a' to accept/lock the full model first.")
            return

        if not hasattr(self, 'continuum_model') or np.all(np.isnan(self.continuum_model)):
             print("No valid continuum model found.")
             return
             
        # Prepare data
        # _cont: wavelength, continuum
        cont_fname = self.filename + "_cont" # Note: filename might have extension
        # Let's handle extension roughly
        base, ext = os.path.splitext(self.filename)
        cont_fname = f"{base}_cont{ext}" if ext else f"{base}_cont.txt"
        norm_fname = f"{base}_norm{ext}" if ext else f"{base}_norm.txt"
        
        # _cont file
        # Filter NaNs? Usually model should be full. If gaps, they are NaNs.
        # Let's save what we have.
        try:
            # We construct full arrays
            # If NaNs exist (e.g. no fit in some regions?), saving might be tricky.
            # User expectation: "entire wavelength range". So hopefully no NaNs.
            
            # _cont
            # columns: wave, cont
            data_cont = np.column_stack((self.wavelength, self.continuum_model))
            np.savetxt(cont_fname, data_cont, header="Wavelength Continuum")
            
            # _norm
            # columns: wave, flux/cont, err/cont
            norm_flux = self.flux / self.continuum_model
            norm_err = self.error / self.continuum_model
            data_norm = np.column_stack((self.wavelength, norm_flux, norm_err))
            np.savetxt(norm_fname, data_norm, header="Wavelength NormFlux NormError")
            
            print(f"Saved fit to:\n  {cont_fname}\n  {norm_fname}")
            
        except Exception as e:
            print(f"Error saving files: {e}")

    def show_help(self, event):
        help_text = """
Spectra Viewer - Help & Controls

Navigation:
  i / o     : Zoom in / out (centered on cursor)
  [ / ]     : Pan left / right
  b / t     : Set bottom / top Y-limit to cursor
  w         : Reset View (Home)
  q         : Quit Application
  l / L     : Toggle Log Scale (l=Y, L=X)

Analysis:
  m         : Identify line at cursor (opens selection)
  z         : Open Line Check Window (3x2 grid of lines)
  Redshift  : Enter value in box to shift lines

Continuum Fitting Mode:
  1. Toggle 'Continuum' checkbox to enable.
  2. Set Polynomial Order.

  Auto Cont : Initialize anchors using sigma-clipped median (100px chunks).
  x         : Add Anchor Point at cursor
  d         : Delete closest Anchor Point
  c         : Compute/Update Polynomial Fit

Equivalent Width:
  e (x2)    : Define integration window.
              - Uses current continuum model if available.
              - Fallback (if no fit):
                - 'Linear Cont' Checked: Linear interpolation.
                - 'Linear Cont' Unchecked: Assumes Cont=1.0.
  
Frozen Segments:
  Space (x2): Freeze a region. Press once to start, again to end.
              - The fit in this region is locked (grey).
              - Anchors inside are excluded from new fits.
              - Overlapping frozen regions are Averaged.
  u         : Unfreeze segment at cursor.
  a         : Accept All. Freezes the current fit in all remaining gaps.
              - Required before Saving.

Saving:
  Save Fit  : Button enabled after 'Accept All'.
              - Saves '_cont.txt' (Continuum)
              - Saves '_norm.txt' (Normalized Spectrum)
"""
        # Create popup window
        hfig = plt.figure(figsize=(6, 9)) # Taller for more text
        hfig.canvas.manager.set_window_title("Spectra Viewer Help")
        hfig.text(0.05, 0.95, help_text, fontsize=10, family='monospace', verticalalignment='top')
        plt.show()

    def handle_ew(self, x):
        if self.ew_start is None:
            # Start
            self.ew_start = x
            self.ew_guide = self.ax.axvline(x, color='green', linestyle='--', alpha=0.5)
            print(f"EW Integration Start: {x:.2f}. Press 'e' again to end.")
            
            # Reset display
            self.text_ew_obs.set_text("Obs: --")
            self.text_ew_rest.set_text("Rest: --")
            if self.ew_poly:
                self.ew_poly.remove()
                self.ew_poly = None
            self.fig.canvas.draw_idle()
        else:
            # End
            start = min(self.ew_start, x)
            end = max(self.ew_start, x)
            print(f"EW Integration Range: {start:.2f} - {end:.2f}")
            
            # Remove guide
            if self.ew_guide:
                self.ew_guide.remove()
                self.ew_guide = None
            self.ew_start = None
            
            # Compute EW
            mask = (self.wavelength >= start) & (self.wavelength <= end)
            if not np.any(mask):
                print("No data in range.")
                return
            
            wave_seg = self.wavelength[mask]
            flux_seg = self.flux[mask]
            
            # Determine Continuum
            cont_seg = None
            
            # 1. Try global model (e.g. from Save/Accept All)
            if hasattr(self, 'continuum_model') and not np.all(np.isnan(self.continuum_model[mask])):
                 cont_seg = self.continuum_model[mask]
                 mode = "Global Model"
                 
            # 2. Try current active poly
            elif hasattr(self, 'current_poly') and self.current_poly is not None:
                 # Check if poly covers this region (it's global usually)
                 cont_seg = self.current_poly(wave_seg)
                 mode = "Active Poly"
            
            # 3. Fallback
            else:
                 if self.linear_cont_fallback:
                     # Interpolate between flux at edges
                     # Find flux closest to start/end
                     idx_s = np.abs(self.wavelength - start).argmin()
                     idx_e = np.abs(self.wavelength - end).argmin()
                     
                     # Simple linear function
                     y1 = self.flux[idx_s]
                     y2 = self.flux[idx_e]
                     x1 = self.wavelength[idx_s]
                     x2 = self.wavelength[idx_e]
                     
                     gradient = (y2 - y1) / (x2 - x1) if x2 != x1 else 0
                     cont_seg = y1 + gradient * (wave_seg - x1)
                     mode = "Linear Fallback"
                 else:
                     # Assume input is already normalized (Cont=1.0)
                     cont_seg = np.ones_like(flux_seg)
                     mode = "Normalized (1.0)"
            
            # Integrate: Integral (1 - F/C) dl
            norm_flux = flux_seg / cont_seg
            integrand = 1.0 - norm_flux
            
            # Trapezoidal integration
            ew_obs = np.trapz(integrand, x=wave_seg)
            self.ew_obs_value = ew_obs # Store for updates
            ew_rest = ew_obs / (1 + self.current_z)
            
            print(f"EW ({mode}): Obs={ew_obs:.3f}, Rest={ew_rest:.3f}")
            
            # Update labels
            self.text_ew_obs.set_text(f"Obs: {ew_obs:.3f} A")
            self.text_ew_rest.set_text(f"Rest: {ew_rest:.3f} A")
            
            # Visualize integration area
            # Fill between continuum and flux
            # Typically show what was integrated
            # We want to shade the area between cont_seg and flux_seg where cont > flux (absorption)
            # But EW technically integrates everything.
            # Let's shade between cont and flux.
            
            self.ew_poly = self.ax.fill_between(wave_seg, cont_seg, flux_seg, color='green', alpha=0.3)
            self.fig.canvas.draw_idle()

    def auto_continuum(self, event):
        if not self.continuum_mode:
            print("Please enable Continuum Mode first.")
            return

        # Clear existing anchors
        # But wait, should we clear only unfrozen anchors?
        # "Initialized anchor points" implies a fresh start usually.
        # But if user has frozen segments, we probably shouldn't mess with them?
        # Prompt: "initialize anchor points".
        # Safe bet: Clear active anchors, respect frozen segments?
        # Or just clear all active anchors.
        
        # Let's clear active anchors (continuum_anchors list).
        for pt in self.continuum_plots:
            pt.remove()
        self.continuum_anchors = []
        self.continuum_plots = []
        
        print("Auto-generating anchors...")
        
        # Chunk size
        chunk_size = 100
        n_pixels = len(self.wavelength)
        
        for i in range(0, n_pixels, chunk_size):
            end = min(i + chunk_size, n_pixels)
            if end - i < 10: continue # Skip tiny chunks
            
            wave_chunk = self.wavelength[i:end]
            flux_chunk = self.flux[i:end]
            
            # Sigma clip stats
            # Ignore NaNs/Infs
            # sigma_clipped_stats returns (mean, median, stddev)
            # We use median.
            mean, median, std = sigma_clipped_stats(flux_chunk, sigma=3.0, maxiters=5)
            
            if np.ma.is_masked(median):
                 # if everything masked, skip
                 continue
                 
            # Add anchor at middle of chunk wavelenght? Or mean wavelength?
            anchor_x = np.mean(wave_chunk)
            anchor_y = float(median) # Ensure float
            
            # Check if this anchor falls into a frozen segment
            # If so, we skip adding it to active anchors
            in_frozen = False
            for seg in self.frozen_segments:
                if seg['range'][0] <= anchor_x <= seg['range'][1]:
                    in_frozen = True
                    break
            
            if not in_frozen:
                self.add_anchor(anchor_x, anchor_y)
                
        print(f"Added {len(self.continuum_anchors)} anchors.")
        self.update_fit() # Optional: Trigger fit immediately? Yes, good UX.

    def add_anchor(self, x, y):
        # Check if in frozen region
        for seg in self.frozen_segments:
            if seg['range'][0] <= x <= seg['range'][1]:
                print("Cannot add anchor in frozen region!")
                return
                
        self.continuum_anchors.append((x, y))
        # Plot pink asterisk
        pt, = self.ax.plot(x, y, 'm*', markersize=10)
        self.continuum_plots.append(pt)
        self.fig.canvas.draw_idle()

    def delete_anchor(self, x, y):
        # Check if in frozen region? Actually user clicks near anchor.
        # Should we prevent deleting if anchor is "under" a frozen region?
        # Yes, prompt says "user cannot add/delete points from frozen regions".
        # But wait, we check cursor position or anchor position?
        # Ideally anchor position.
        
        if not self.continuum_anchors: return
        
        # Find closest
        dists = [np.sqrt((ax - x)**2 + (ay - y)**2) for ax, ay in self.continuum_anchors]
        idx = np.argmin(dists)
        anchor = self.continuum_anchors[idx]
        
        # Check if anchor is in frozen region
        for seg in self.frozen_segments:
            if seg['range'][0] <= anchor[0] <= seg['range'][1]:
                print("Cannot delete anchor in frozen region!")
                return
        
        # Remove from list and plot
        self.continuum_anchors.pop(idx)
        pt = self.continuum_plots.pop(idx)
        pt.remove()
        
        self.fig.canvas.draw_idle()

    def update_fit(self):
        # Exclude anchors in frozen regions
        active_anchors = []
        for p in self.continuum_anchors:
            in_frozen = False
            for seg in self.frozen_segments:
                if seg['range'][0] <= p[0] <= seg['range'][1]:
                    in_frozen = True
                    break
            if not in_frozen:
                active_anchors.append(p)
                
        if len(active_anchors) < self.poly_order + 1:
            print(f"Not enough active anchors for order {self.poly_order} fit.")
            return

        # Sort anchors by x
        anchors = sorted(active_anchors, key=lambda p: p[0])
        x_anchors = np.array([p[0] for p in anchors])
        y_anchors = np.array([p[1] for p in anchors])
        
        try:
            # Polyfit
            coeffs = np.polyfit(x_anchors, y_anchors, self.poly_order)
            self.current_poly = np.poly1d(coeffs) # Store Poly1d obj
            
            # Generate line on view grid for display? Or full grid?
            # Let's compute on full wavelength grid but only plot within view or handle appropriately?
            # Plotting usually uses view limits for speed, but for consistency let's just plot on x_plot within view.
            
            xmin, xmax = self.ax.get_xlim()
            x_plot = np.linspace(xmin, xmax, 1000)
            y_plot = self.current_poly(x_plot)
            
            # Remove old fit
            if self.continuum_fit_line:
                self.continuum_fit_line.remove()
                
            self.continuum_fit_line, = self.ax.plot(x_plot, y_plot, 'k:', linewidth=2)
            self.fig.canvas.draw_idle()
            
        except Exception as e:
            print(f"Fit failed: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python spectra_viewer.py <filename>")
        sys.exit(1)
    
    filename = sys.argv[1]
    viewer = SpectraViewer(filename)
