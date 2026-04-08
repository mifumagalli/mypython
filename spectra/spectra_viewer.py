import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import UnivariateSpline


def load_spectrum(filename):
    """Loads spectrum from an ASCII file (wavelength, flux, error)."""
    try:
        data = np.loadtxt(filename)
        if data.shape[1] < 2:
            raise ValueError("File must have at least 2 columns (wavelength, flux)")
        wavelength = data[:, 0]
        flux = data[:, 1]
        error = data[:, 2] if data.shape[1] > 2 else np.zeros_like(flux)
        mask = np.logical_not(np.isfinite(flux))
        flux[mask] = 0
        error[mask] = 0
        return wavelength, flux, error
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)


from matplotlib.widgets import Button, TextBox, CheckButtons, RadioButtons
from linetools.lists.linelist import LineList
from astropy.stats import sigma_clipped_stats
import warnings

warnings.filterwarnings("ignore")


class SpectraViewer:
    def __init__(self, filename):
        self.filename = filename
        self.wavelength, self.flux, self.error = load_spectrum(filename)
        self.linelist = LineList("Strong")
        self.lines_visible = False
        self.current_z = 0.0
        self.vlines = []

        self.line_ew_state = {}
        self.fig = plt.figure(figsize=(12, 6))
        self.gs = self.fig.add_gridspec(
            2, 2, width_ratios=[4, 1], height_ratios=[4, 1], hspace=0.1
        )

        self.ax = self.fig.add_subplot(self.gs[:, 0])
        (self.line,) = self.ax.plot(self.wavelength, self.flux, "b-", label="Flux")
        (self.error_line,) = self.ax.plot(
            self.wavelength, self.error, "orange", label="Error", alpha=0.7
        )

        self.ax.set_xlabel("Wavelength")
        self.ax.set_ylabel("Flux")
        self.ax.set_title(f"Spectrum: {filename}")
        self.ax.legend()

        self.ax.autoscale(enable=True, axis="both", tight=True)
        self.ax.relim()
        self.ax.autoscale_view()

        self.panel_ax = self.fig.add_subplot(self.gs[:, 1])
        self.panel_ax.axis("off")

        # --- Spectrum Controls ---
        self.fig.text(0.82, 0.94, "Spectrum Controls", fontweight="bold", fontsize=10)
        line1 = plt.Line2D(
            [0.82, 0.97],
            [0.92, 0.92],
            transform=self.fig.transFigure,
            color="k",
            linewidth=1,
        )
        self.fig.add_artist(line1)

        ax_check = self.fig.add_axes([0.82, 0.85, 0.15, 0.05])
        self.check = CheckButtons(ax_check, ["Show Lines"], [self.lines_visible])
        self.check.on_clicked(self.toggle_lines)

        ax_text = self.fig.add_axes([0.875, 0.75, 0.095, 0.05])
        self.text_box = TextBox(ax_text, "Redshift: ", initial=str(self.current_z))
        self.text_box.on_submit(self.submit_redshift)

        # --- Continuum Fit ---
        self.fig.text(0.82, 0.66, "Continuum Fit", fontweight="bold", fontsize=10)
        line2 = plt.Line2D(
            [0.82, 0.97],
            [0.64, 0.64],
            transform=self.fig.transFigure,
            color="k",
            linewidth=1,
        )
        self.fig.add_artist(line2)

        self.continuum_mode = False
        self.continuum_anchors = []
        self.continuum_plots = []
        self.continuum_fit_line = None
        self.poly_order = 1
        self.current_poly = None

        self.frozen_segments = []
        self.freeze_start = None
        self.freeze_guide = None
        self.model_accepted = False
        self.frozen_line_obj = None

        ax_check_cont = self.fig.add_axes([0.82, 0.57, 0.15, 0.05])
        self.check_cont = CheckButtons(
            ax_check_cont, ["Continuum"], [self.continuum_mode]
        )
        self.check_cont.on_clicked(self.toggle_continuum_mode)

        ax_text_poly = self.fig.add_axes([0.89, 0.49, 0.08, 0.05])
        self.text_poly = TextBox(
            ax_text_poly, "Poly Order: ", initial=str(self.poly_order)
        )
        self.text_poly.on_submit(self.submit_poly_order)

        # Make Anchors (from flux sigma-clip)
        ax_auto = self.fig.add_axes([0.82, 0.42, 0.15, 0.05])
        self.btn_auto = Button(ax_auto, "Make Anchors")
        self.btn_auto.on_clicked(self.auto_continuum)

        # Fit Spline
        ax_spline = self.fig.add_axes([0.82, 0.36, 0.15, 0.05])
        self.btn_spline = Button(ax_spline, "Fit Spline")
        self.btn_spline.on_clicked(self.fit_spline)

        # Load Continuum (file dialog)
        #ax_load_cont = self.fig.add_axes([0.82, 0.30, 0.15, 0.05])
        #self.btn_load_cont = Button(ax_load_cont, "Load Continuum")
        #self.btn_load_cont.on_clicked(self.load_continuum)
        #self.btn_load_cont.label.set_fontsize(8)

        # -----------------------------------------------------------------
        # Auto Continuum button (NEW)
        # Finds <basename>_cont<ext> automatically, loads it, places anchors
        # on the continuum values, enables continuum mode, refreshes ratio plot.
        # Button is tinted green when the companion file already exists.
        # -----------------------------------------------------------------
        ax_auto_cont = self.fig.add_axes([0.82, 0.30, 0.15, 0.05])
        self.btn_auto_cont = Button(ax_auto_cont, "Load Continuum")
        self.btn_auto_cont.on_clicked(self.auto_load_and_anchor_continuum)
        self.btn_auto_cont.label.set_fontsize(8)
        self._update_auto_cont_button_color()
        # -----------------------------------------------------------------

        # Accept / Save
        ax_save = self.fig.add_axes([0.82, 0.24, 0.15, 0.05])
        self.btn_save = Button(ax_save, "Accept fit to save")
        self.btn_save.on_clicked(self.save_fit)
        self.btn_save.label.set_fontsize(8)
        self.btn_save.color = "0.9"
        self.btn_save.hovercolor = "0.9"

        # --- Equivalent Width ---
        self.fig.text(0.82, 0.21, "Equivalent Width", fontweight="bold", fontsize=10)
        line3 = plt.Line2D(
            [0.82, 0.97],
            [0.195, 0.195],
            transform=self.fig.transFigure,
            color="k",
            linewidth=1,
        )
        self.fig.add_artist(line3)

        ax_check_ew = self.fig.add_axes([0.82, 0.15, 0.15, 0.03])
        self.linear_cont_fallback = False
        self.check_ew = CheckButtons(
            ax_check_ew, ["Linear Cont"], [self.linear_cont_fallback]
        )
        self.check_ew.on_clicked(self.toggle_linear_fallback)

        self.text_ew_obs = self.fig.text(0.82, 0.13, "Obs: --", fontsize=9)
        self.text_ew_rest = self.fig.text(0.82, 0.11, "Rest: --", fontsize=9)

        self.ew_start = None
        self.ew_guide = None
        self.ew_poly = None
        self.ew_obs_value = None

        ax_help = self.fig.add_axes([0.82, 0.02, 0.07, 0.025])
        self.btn_help = Button(ax_help, "Help")
        self.btn_help.label.set_fontsize(7)
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

    # ------------------------------------------------------------------
    # Helpers for Auto Continuum
    # ------------------------------------------------------------------

    def _cont_filename(self):
        """Return the expected companion _cont filename."""
        base, ext = os.path.splitext(self.filename)
        return f"{base}_cont{ext}" if ext else f"{base}_cont.txt"

    def _update_auto_cont_button_color(self):
        """Green tint when companion file exists, grey otherwise."""
        if os.path.exists(self._cont_filename()):
            self.btn_auto_cont.color = "palegreen"
            self.btn_auto_cont.hovercolor = "lightgreen"
        else:
            self.btn_auto_cont.color = "0.85"
            self.btn_auto_cont.hovercolor = "0.95"

    # ------------------------------------------------------------------
    # Auto Continuum: load _cont file → anchors → ratio plot  (NEW)
    # ------------------------------------------------------------------

    def auto_load_and_anchor_continuum(self, event):
        """
        1. Locate <basename>_cont<ext> next to the spectrum file.
        2. Load + interpolate it onto the spectrum wavelength grid.
        3. Display it as a grey continuum line on the main axes.
        4. Clear any old anchors and place new ones (one per 100-px chunk)
           sampled from the continuum values (not the raw flux).
        5. Enable continuum mode if it isn't already on.
        6. Refresh the ratio plot.
        """
        cont_file = self._cont_filename()

        if not os.path.exists(cont_file):
            print(
                f"Auto Continuum: companion file not found.\n"
                f"  Expected: {cont_file}\n"
                f"  Use 'Load Continuum' to choose a file manually."
            )
            return

        # 1. Load -------------------------------------------------------
        try:
            data = np.loadtxt(cont_file, comments="#")
        except Exception as e:
            print(f"Auto Continuum: could not read '{cont_file}': {e}")
            return

        if data.ndim != 2 or data.shape[1] < 2:
            print(
                "Auto Continuum: file must have at least 2 columns "
                "(wavelength, continuum)."
            )
            return

        cont_wave = data[:, 0]
        cont_flux = data[:, 1]

        sort_idx = np.argsort(cont_wave)
        cont_wave = cont_wave[sort_idx]
        cont_flux = cont_flux[sort_idx]

        spec_min, spec_max = self.wavelength.min(), self.wavelength.max()
        cont_min, cont_max = cont_wave.min(), cont_wave.max()

        if cont_max < spec_min or cont_min > spec_max:
            print(
                f"Auto Continuum: wavelength ranges do not overlap.\n"
                f"  Continuum: [{cont_min:.1f}, {cont_max:.1f}]\n"
                f"  Spectrum:  [{spec_min:.1f}, {spec_max:.1f}]"
            )
            return

        # 2. Interpolate onto spectrum grid ------------------------------
        self.continuum_model = np.interp(
            self.wavelength, cont_wave, cont_flux, left=np.nan, right=np.nan
        )

        # 3. Display continuum line -------------------------------------
        if self.frozen_line_obj is not None:
            try:
                self.frozen_line_obj.remove()
            except Exception:
                pass

        (self.frozen_line_obj,) = self.ax.plot(
            self.wavelength,
            self.continuum_model,
            color="grey",
            linewidth=2,
            linestyle="-",
            label="Auto Continuum",
        )

        # Clear any active spline so the loaded model takes precedence
        if self.continuum_fit_line is not None:
            try:
                self.continuum_fit_line.remove()
            except Exception:
                pass
            self.continuum_fit_line = None
        self.current_poly = None

        # Mark model as accepted → Save button becomes active
        self.model_accepted = True
        self.btn_save.label.set_text("Save Fit")
        self.btn_save.color = "0.85"
        self.btn_save.hovercolor = "0.95"

        # 4. Generate anchors from the continuum ------------------------
        for pt in self.continuum_plots:
            try:
                pt.remove()
            except Exception:
                pass
        self.continuum_anchors = []
        self.continuum_plots = []

        anchor_count = 0
        has_anchor_col = False
        if data.shape[1] >= 3:
            is_anchor = data[:, -1]
            if set(np.unique(is_anchor)).issubset({0, 1}):
                is_anchor_sorted = is_anchor[sort_idx]
                has_anchor_col = True

        if has_anchor_col:
            for wave, cont, is_anch in zip(cont_wave, cont_flux, is_anchor_sorted):
                if is_anch == 1 and np.isfinite(cont):
                    in_frozen = any(
                        seg["range"][0] <= wave <= seg["range"][1]
                        for seg in self.frozen_segments
                    )
                    if not in_frozen:
                        self.continuum_anchors.append((wave, cont))
                        (pt,) = self.ax.plot(wave, cont, "m*", markersize=10)
                        self.continuum_plots.append(pt)
                        anchor_count += 1
        else:
            chunk_size = 100
            n_pixels = len(self.wavelength)

            for i in range(0, n_pixels, chunk_size):
                end = min(i + chunk_size, n_pixels)

                wave_chunk = self.wavelength[i:end]
                cont_chunk = self.continuum_model[i:end]

                finite_mask = np.isfinite(cont_chunk)
                if np.sum(finite_mask) < 5:
                    continue

                anchor_x = float(np.mean(wave_chunk[finite_mask]))
                anchor_y = float(np.mean(cont_chunk[finite_mask]))

                in_frozen = any(
                    seg["range"][0] <= anchor_x <= seg["range"][1]
                    for seg in self.frozen_segments
                )
                if in_frozen:
                    continue

                self.continuum_anchors.append((anchor_x, anchor_y))
                (pt,) = self.ax.plot(anchor_x, anchor_y, "m*", markersize=10)
                self.continuum_plots.append(pt)
                anchor_count += 1

        n_valid = int(np.sum(np.isfinite(self.continuum_model)))
        print(
            f"Auto Continuum: loaded '{os.path.basename(cont_file)}' "
            f"({n_valid}/{len(self.wavelength)} px covered), "
            f"placed {anchor_count} anchors"
            + (" (from IsAnchor column)." if has_anchor_col else " (auto-generated).")
        )

        # 5. Enable continuum mode if not already on --------------------
        if not self.continuum_mode:
            self.check_cont.set_active(0)  # triggers toggle_continuum_mode
            # toggle_continuum_mode calls update_ratio_plot, so we are done.
        else:
            # Mode already on — just refresh the ratio panel.
            self.update_ratio_plot()

        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------------

    def connect(self):
        self.cid = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.ax.callbacks.connect("xlim_changed", self.on_xlim_changed)

    def on_xlim_changed(self, event_ax):
        if self.lines_visible:
            self.draw_lines()

    def toggle_lines(self, label):
        self.lines_visible = not self.lines_visible
        self.draw_lines()

    def update_ratio_plot(self):
        if not hasattr(self, "ax_ratio") or not self.ax_ratio.get_visible():
            return

        current_xlim = self.ax.get_xlim()
        self.ax_ratio.clear()
        self.ax_ratio.axhline(1, color="k", linestyle="-")

        # Build best available continuum:
        # Priority 1 — loaded continuum_model (from auto or load buttons)
        # Priority 2 — active spline/poly
        # Priority 3 — frozen segments override both where they exist
        cont = np.full_like(self.flux, np.nan)

        if hasattr(self, "continuum_model") and not np.all(
            np.isnan(self.continuum_model)
        ):
            cont = self.continuum_model.copy()
        elif self.current_poly is not None:
            cont = self.current_poly(self.wavelength)

        if self.frozen_segments:
            sum_y = np.zeros(len(self.wavelength))
            count_y = np.zeros(len(self.wavelength))
            for seg in self.frozen_segments:
                idx = seg["mask_idx"]
                sum_y[idx] += seg["y"]
                count_y[idx] += 1
            mask_valid = count_y > 0
            cont[mask_valid] = sum_y[mask_valid] / count_y[mask_valid]

        if not np.all(np.isnan(cont)):
            ratio = self.flux / cont
            ratio_err = self.error / cont
            self.ax_ratio.plot(
                self.wavelength, ratio, "b-", label="Flux / Cont", linewidth=1
            )
            self.ax_ratio.fill_between(
                self.wavelength,
                ratio - ratio_err,
                ratio + ratio_err,
                color="orange",
                alpha=0.5,
                label="+/- 1 Sigma",
            )
            ratio_finite = ratio[np.isfinite(ratio)]
            if len(ratio_finite) > 0:
                _, _, stdVal = sigma_clipped_stats(ratio_finite, sigma=3.0, maxiters=5)
                stdVal = max(stdVal, 0.01)
                self.ax_ratio.set_ylim(1 - 3 * stdVal, 1 + 3 * stdVal)

        self.ax_ratio.set_ylabel("Flux/Cont")
        self.ax_ratio.set_xlabel("Wavelength")
        self.ax_ratio.set_xlim(current_xlim)
        self.fig.canvas.draw_idle()

    def toggle_continuum_mode(self, label):
        self.continuum_mode = not self.continuum_mode
        if self.continuum_mode:
            self.ax.set_subplotspec(self.gs[0, 0])
            self.ax.set_position(self.gs[0, 0].get_position(self.fig))
            self.ax.tick_params(labelbottom=False)
            self.ax.set_xlabel("")

            if not hasattr(self, "ax_ratio"):
                self.ax_ratio = self.fig.add_subplot(self.gs[1, 0], sharex=self.ax)
            else:
                self.ax_ratio.set_visible(True)

            self.update_ratio_plot()
        else:
            self.ax.set_subplotspec(self.gs[:, 0])
            self.ax.set_position(self.gs[:, 0].get_position(self.fig))
            self.ax.tick_params(labelbottom=True)
            self.ax.set_xlabel("Wavelength")
            if hasattr(self, "ax_ratio"):
                self.ax_ratio.set_visible(False)

        self.fig.canvas.draw_idle()

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
            if self.ew_obs_value is not None:
                ew_rest = self.ew_obs_value / (1 + self.current_z)
                self.text_ew_rest.set_text(f"Rest: {ew_rest:.3f} A")
                self.fig.canvas.draw_idle()
        except ValueError:
            print("Invalid redshift value")

    def draw_lines(self):
        for line in self.vlines:
            line.remove()
        self.vlines = []

        if not self.lines_visible:
            self.fig.canvas.draw_idle()
            return

        xlim = self.ax.get_xlim()
        observed_waves = self.linelist.wrest.value * (1 + self.current_z)
        mask = (observed_waves > xlim[0]) & (observed_waves < xlim[1])
        visible_waves = observed_waves[mask]
        visible_names = self.linelist.name[mask]

        for wv, name in zip(visible_waves, visible_names):
            vline = self.ax.axvline(wv, color="k", linestyle="--", alpha=0.5)
            self.vlines.append(vline)
            text = self.ax.text(
                wv,
                0.95,
                name,
                rotation=90,
                verticalalignment="top",
                horizontalalignment="right",
                transform=self.ax.get_xaxis_transform(),
                fontsize=8,
                clip_on=True,
            )
            self.vlines.append(text)

        self.fig.canvas.draw_idle()

    def identify_line(self, cursor_x):
        sorted_indices = self.linelist.wrest.argsort()
        sorted_names = self.linelist.name[sorted_indices]
        sorted_wrest = self.linelist.wrest[sorted_indices]

        labels = [
            f"{name} ({wrest.value:.1f})"
            for name, wrest in zip(sorted_names, sorted_wrest)
        ]

        N_per_col = 15
        COLS_PER_PAGE = 3
        width = 2.5 * COLS_PER_PAGE
        height = min(6, 0.4 * N_per_col) + 1.0

        id_fig = plt.figure(figsize=(width, height))
        id_fig.canvas.manager.set_window_title("Select Line")

        self.current_page = 0
        self.radio_buttons = []
        self.page_buttons = []

        def draw_page():
            id_fig.clf()
            self.radio_buttons = []
            self.page_buttons = []

            start_x = 0.05
            col_width = 0.9 / COLS_PER_PAGE
            page_start_idx = self.current_page * N_per_col * COLS_PER_PAGE

            for i in range(COLS_PER_PAGE):
                col_offset = i * N_per_col
                start = page_start_idx + col_offset
                if start >= len(labels):
                    break
                end = min(start + N_per_col, len(labels))
                chunk_labels = labels[start:end]

                ax_col = id_fig.add_axes(
                    [start_x + i * col_width, 0.15, col_width, 0.8]
                )
                ax_col.axis("off")
                radio = RadioButtons(ax_col, chunk_labels, active=0)
                self.radio_buttons.append(radio)

                def make_callback(offset):
                    def on_select(label):
                        if label not in labels:
                            return
                        idx = labels.index(label)
                        rest_wave = sorted_wrest[idx].value
                        new_z = (cursor_x / rest_wave) - 1
                        print(
                            f"Selected {sorted_names[idx]} at rest {rest_wave:.2f}. "
                            f"Obs: {cursor_x:.2f}. New z: {new_z:.4f}"
                        )
                        self.current_z = new_z
                        self.text_box.set_val(f"{new_z:.4f}")
                        if not self.lines_visible:
                            self.check.set_active(0)
                        plt.close(id_fig)

                    return on_select

                radio.on_clicked(make_callback(start))

            ax_prev = id_fig.add_axes([0.1, 0.02, 0.3, 0.08])
            b_prev = Button(ax_prev, "Previous")
            ax_next = id_fig.add_axes([0.6, 0.02, 0.3, 0.08])
            b_next = Button(ax_next, "Next")
            self.page_buttons.extend([b_prev, b_next])

            def prev_clicked(event):
                if self.current_page > 0:
                    self.current_page -= 1
                    draw_page()

            def next_clicked(event):
                if (self.current_page + 1) * N_per_col * COLS_PER_PAGE < len(labels):
                    self.current_page += 1
                    draw_page()

            b_prev.on_clicked(prev_clicked)
            b_next.on_clicked(next_clicked)
            id_fig.canvas.draw()

        draw_page()
        plt.show()

    def show_line_check_window(self):
        min_obs = np.min(self.wavelength)
        max_obs = np.max(self.wavelength)
        all_obs_waves = self.linelist.wrest.value * (1 + self.current_z)
        in_range_mask = (all_obs_waves >= min_obs) & (all_obs_waves <= max_obs)

        if not np.any(in_range_mask):
            print(
                f"No lines in range {min_obs:.1f} - {max_obs:.1f} A for z={self.current_z:.4f}"
            )
            return

        filtered_names = self.linelist.name[in_range_mask]
        filtered_wrest = self.linelist.wrest[in_range_mask]
        sorted_indices = filtered_wrest.argsort()
        sorted_names = filtered_names[sorted_indices]
        sorted_wrest = filtered_wrest[sorted_indices]

        COLS, ROWS = 3, 2
        N_PER_PAGE = COLS * ROWS

        check_fig = plt.figure(figsize=(15, 8))
        check_fig.canvas.manager.set_window_title(
            f"Line Check (z={self.current_z:.4f})"
        )

        self.check_current_page = 0
        self.check_buttons = []

        def draw_check_page():
            check_fig.clf()
            self.check_buttons = []
            self.check_axes_map = {}

            gs = check_fig.add_gridspec(
                ROWS, COLS, bottom=0.15, top=0.95, hspace=0.4, wspace=0.3
            )
            start_idx = self.check_current_page * N_PER_PAGE

            for i in range(N_PER_PAGE):
                idx = start_idx + i
                if idx >= len(sorted_names):
                    break

                rest_wave = sorted_wrest[idx].value
                name = sorted_names[idx]
                obs_wave = rest_wave * (1 + self.current_z)

                c_kms = 299792.458
                delta_wave = obs_wave * (1000.0 / c_kms)
                w_min = obs_wave - delta_wave
                w_max = obs_wave + delta_wave

                mask = (self.wavelength >= w_min) & (self.wavelength <= w_max)
                row = i // COLS
                col = i % COLS
                ax = check_fig.add_subplot(gs[row, col])
                self.check_axes_map[ax] = (name, rest_wave)

                if np.any(mask):
                    ax.plot(self.wavelength[mask], self.flux[mask], "b-", linewidth=1)
                    ax.plot(
                        self.wavelength[mask],
                        self.error[mask],
                        "orange",
                        alpha=0.5,
                        linewidth=1,
                    )

                ax.axvline(obs_wave, color="k", linestyle="--")
                ax.set_title(
                    f"{name}\nRest: {rest_wave:.1f}, Obs: {obs_wave:.1f}", fontsize=10
                )
                ax.set_xlim(w_min, w_max)

                if np.any(mask):
                    y_seg = self.flux[mask]
                    if len(y_seg) > 0:
                        ymin, ymax = np.min(y_seg), np.max(y_seg)
                        margin = (ymax - ymin) * 0.1 if ymax != ymin else 1.0
                        ax.set_ylim(ymin - margin, ymax + margin)

            if self.check_current_page > 0:
                ax_prev = check_fig.add_axes([0.1, 0.05, 0.1, 0.05])
                b_prev = Button(ax_prev, "Previous")
                b_prev.on_clicked(lambda e: change_page(-1))
                self.check_buttons.append(b_prev)

            if (self.check_current_page + 1) * N_PER_PAGE < len(sorted_names):
                ax_next = check_fig.add_axes([0.8, 0.05, 0.1, 0.05])
                b_next = Button(ax_next, "Next")
                b_next.on_clicked(lambda e: change_page(1))
                self.check_buttons.append(b_next)

            ax_measure = check_fig.add_axes([0.65, 0.05, 0.12, 0.05])
            b_measure = Button(ax_measure, "Save EW table")
            b_measure.on_clicked(self.measure_ew_batch)
            self.check_buttons.append(b_measure)

            check_fig.text(
                0.5,
                0.075,
                f"Page {self.check_current_page + 1} / {int(np.ceil(len(sorted_names)/N_PER_PAGE))}",
                horizontalalignment="center",
                fontsize=12,
            )
            check_fig.text(
                0.5,
                0.02,
                "Controls: 'u' Select | 'd' Deselect | 'e'+'e' Custom Window",
                horizontalalignment="center",
                fontsize=10,
                color="blue",
            )
            check_fig.canvas.draw()

        def change_page(delta):
            self.check_current_page += delta
            draw_check_page()

        check_fig.canvas.mpl_connect("key_press_event", self.on_check_key)
        self.check_fig = check_fig
        self.sorted_names_cache = sorted_names
        self.sorted_wrest_cache = sorted_wrest
        self.check_axes_cache = {}
        draw_check_page()
        plt.show()

    def on_check_key(self, event):
        if not event.inaxes:
            return
        if (
            not hasattr(self, "check_axes_map")
            or event.inaxes not in self.check_axes_map
        ):
            return

        name, wrest_val = self.check_axes_map[event.inaxes]
        key_id = (name, wrest_val)

        if key_id not in self.line_ew_state:
            self.line_ew_state[key_id] = {"selected": False, "custom_window": None}

        state = self.line_ew_state[key_id]

        if event.key == "u":
            state["selected"] = True
            self.update_check_panel_title(event.inaxes, name, wrest_val, state)
        elif event.key == "d":
            state["selected"] = False
            state["custom_window"] = None
            self.update_check_panel_title(event.inaxes, name, wrest_val, state)
        elif event.key == "e":
            if "ew_wait_start" not in state:
                state["ew_wait_start"] = event.xdata
                guide = event.inaxes.axvline(event.xdata, color="green", linestyle="--")
                state["guide_artist"] = guide
                event.inaxes.set_title("Define End...", color="blue")
                event.canvas.draw_idle()
            else:
                start = min(state["ew_wait_start"], event.xdata)
                end = max(state["ew_wait_start"], event.xdata)
                state["custom_window"] = (start, end)
                if "guide_artist" in state:
                    state["guide_artist"].remove()
                    del state["guide_artist"]
                del state["ew_wait_start"]
                state["selected"] = True
                self.draw_check_panel_extras(event.inaxes, state)
                self.update_check_panel_title(event.inaxes, name, wrest_val, state)

    def update_check_panel_title(self, ax, name, wrest, state):
        obs = wrest * (1 + self.current_z)
        title_text = f"{name}\nRest: {wrest:.1f}, Obs: {obs:.1f}"
        color = "black"
        if state.get("selected"):
            color = "green"
        if state.get("selected") is False:
            color = "red"
        ax.set_title(title_text, color=color, fontsize=10)
        self.check_fig.canvas.draw_idle()

    def draw_check_panel_extras(self, ax, state):
        if "window_patch" in state and state["window_patch"]:
            try:
                state["window_patch"].remove()
            except Exception:
                pass
            state["window_patch"] = None

        if state.get("custom_window"):
            w_min, w_max = state["custom_window"]
            patch = ax.axvspan(w_min, w_max, color="green", alpha=0.2)
            state["window_patch"] = patch

        self.check_fig.canvas.draw_idle()

    def measure_ew_batch(self, event):
        print("Measuring EW for selected lines...")
        results = []

        for key_id, state in self.line_ew_state.items():
            if state.get("selected"):
                name, wrest = key_id

                if state.get("custom_window"):
                    w_start, w_end = state["custom_window"]
                else:
                    obs_wave = wrest * (1 + self.current_z)
                    c_kms = 299792.458
                    delta = obs_wave * (200.0 / c_kms)
                    w_start = obs_wave - delta
                    w_end = obs_wave + delta

                mask = (self.wavelength >= w_start) & (self.wavelength <= w_end)
                if not np.any(mask):
                    print(f"Skipping {name}: No data in window.")
                    continue

                wave_seg = self.wavelength[mask]
                flux_seg = self.flux[mask]
                err_seg = self.error[mask]

                if hasattr(self, "continuum_model") and not np.all(
                    np.isnan(self.continuum_model[mask])
                ):
                    cont_seg = self.continuum_model[mask]
                elif self.current_poly is not None:
                    cont_seg = self.current_poly(wave_seg)
                else:
                    cont_seg = np.ones_like(flux_seg)

                ew_obs = np.trapz(1.0 - flux_seg / cont_seg, x=wave_seg)
                dl = np.gradient(wave_seg)
                ew_err_obs = np.sqrt(np.sum(((err_seg / cont_seg) * dl) ** 2))
                ew_rest = ew_obs / (1 + self.current_z)
                ew_err_rest = ew_err_obs / (1 + self.current_z)
                obs_center = wrest * (1 + self.current_z)

                results.append(
                    (
                        name,
                        wrest,
                        obs_center,
                        w_start,
                        w_end,
                        ew_obs,
                        ew_err_obs,
                        ew_rest,
                        ew_err_rest,
                    )
                )

        if not results:
            print("No lines selected or valid for measurement.")
            return

        fname = os.path.splitext(self.filename)[0] + "_z{:0.4f}_EW.txt".format(
            self.current_z
        )
        try:
            with open(fname, "w") as f:
                f.write(
                    "# Name RestWave ObsWave WinStart WinEnd EW_Obs Err_Obs EW_Rest Err_Rest\n"
                )
                for res in results:
                    f.write(
                        "{:s} {:.2f} {:.2f} {:.2f} {:.2f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                            *res
                        )
                    )
            print(f"Saved {len(results)} measurements to {fname}")
        except Exception as e:
            print(f"Error saving EW file: {e}")

    # ------------------------------------------------------------------
    # Load Continuum (file dialog)
    # ------------------------------------------------------------------

    def load_continuum(self, event):
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)

        cont_file = filedialog.askopenfilename(
            title="Select continuum file",
            filetypes=[
                ("Text / ASCII files", "*.txt *.dat *.ascii *.fit *.fits"),
                ("All files", "*.*"),
            ],
        )
        root.destroy()

        if not cont_file:
            print("Load Continuum: cancelled.")
            return

        try:
            data = np.loadtxt(cont_file, comments="#")
        except Exception as e:
            print(f"Load Continuum: could not read '{cont_file}': {e}")
            return

        if data.ndim != 2 or data.shape[1] < 2:
            print("Load Continuum: file must have at least 2 columns.")
            return

        cont_wave = data[:, 0]
        cont_flux = data[:, 1]
        sort_idx = np.argsort(cont_wave)
        cont_wave = cont_wave[sort_idx]
        cont_flux = cont_flux[sort_idx]

        spec_min, spec_max = self.wavelength.min(), self.wavelength.max()
        if cont_wave.max() < spec_min or cont_wave.min() > spec_max:
            print("Load Continuum: wavelength ranges do not overlap.")
            return

        self.continuum_model = np.interp(
            self.wavelength, cont_wave, cont_flux, left=np.nan, right=np.nan
        )

        if self.frozen_line_obj is not None:
            try:
                self.frozen_line_obj.remove()
            except Exception:
                pass

        (self.frozen_line_obj,) = self.ax.plot(
            self.wavelength,
            self.continuum_model,
            color="grey",
            linewidth=2,
            linestyle="-",
            label="Loaded Continuum",
        )

        self.model_accepted = True
        self.btn_save.label.set_text("Save Fit")
        self.btn_save.color = "0.85"
        self.btn_save.hovercolor = "0.95"

        # Check for anchors in the loaded file
        for pt in self.continuum_plots:
            try:
                pt.remove()
            except Exception:
                pass
        self.continuum_anchors = []
        self.continuum_plots = []
        
        anchor_count = 0
        has_anchor_col = False
        if data.shape[1] >= 3:
            is_anchor = data[:, -1]
            if set(np.unique(is_anchor)).issubset({0, 1}):
                is_anchor_sorted = is_anchor[sort_idx]
                has_anchor_col = True
                
        if has_anchor_col:
            for wave, cont, is_anch in zip(cont_wave, cont_flux, is_anchor_sorted):
                if is_anch == 1 and np.isfinite(cont):
                    in_frozen = any(
                        seg["range"][0] <= wave <= seg["range"][1]
                        for seg in self.frozen_segments
                    )
                    if not in_frozen:
                        self.continuum_anchors.append((wave, cont))
                        (pt,) = self.ax.plot(wave, cont, "m*", markersize=10)
                        self.continuum_plots.append(pt)
                        anchor_count += 1

        if self.continuum_fit_line is not None:
            try:
                self.continuum_fit_line.remove()
            except Exception:
                pass
            self.continuum_fit_line = None
        self.current_poly = None

        self.update_ratio_plot()
        self.fig.canvas.draw_idle()

        n_valid = int(np.sum(np.isfinite(self.continuum_model)))
        anchor_msg = f", placed {anchor_count} anchors (from IsAnchor column)" if has_anchor_col else ""
        print(
            f"Load Continuum: loaded '{os.path.basename(cont_file)}' — "
            f"{n_valid}/{len(self.wavelength)} pixels covered{anchor_msg}."
        )

    # ------------------------------------------------------------------

    def on_key(self, event):
        if not event.inaxes and event.key != "l":
            if event.key != "l":
                return

        valid_axes = [self.ax]
        if hasattr(self, "ax_ratio") and self.ax_ratio.get_visible():
            valid_axes.append(self.ax_ratio)

        nav_keys = [
            "i",
            "o",
            "[",
            "]",
            "b",
            "t",
            "w",
            "l",
            "L",
            "m",
            "z",
            "e",
            "x",
            "d",
            "c",
            " ",
            "u",
            "a",
        ]
        if event.inaxes not in valid_axes and event.key in nav_keys:
            return

        xlim = self.ax.get_xlim()
        x_range = xlim[1] - xlim[0]
        key = event.key

        if key == "i":
            new_range = x_range * 0.5
            center = event.xdata
            self.ax.set_xlim(center - new_range / 2, center + new_range / 2)
            if hasattr(self, "ax_ratio") and self.ax_ratio.get_visible():
                self.ax_ratio.set_xlim(center - new_range / 2, center + new_range / 2)
        elif key == "o":
            new_range = x_range * 2.0
            center = event.xdata
            self.ax.set_xlim(center - new_range / 2, center + new_range / 2)
            if hasattr(self, "ax_ratio") and self.ax_ratio.get_visible():
                self.ax_ratio.set_xlim(center - new_range / 2, center + new_range / 2)
        elif key == "[":
            pan = x_range * 0.2
            self.ax.set_xlim(xlim[0] - pan, xlim[1] - pan)
            if hasattr(self, "ax_ratio") and self.ax_ratio.get_visible():
                self.ax_ratio.set_xlim(xlim[0] - pan, xlim[1] - pan)
        elif key == "]":
            pan = x_range * 0.2
            self.ax.set_xlim(xlim[0] + pan, xlim[1] + pan)
            if hasattr(self, "ax_ratio") and self.ax_ratio.get_visible():
                self.ax_ratio.set_xlim(xlim[0] + pan, xlim[1] + pan)
        elif key == "b":
            event.inaxes.set_ylim(bottom=event.ydata)
        elif key == "t":
            event.inaxes.set_ylim(top=event.ydata)
        elif key == "l":
            event.inaxes.set_yscale(
                "log" if event.inaxes.get_yscale() == "linear" else "linear"
            )
        elif key == "L":
            event.inaxes.set_xscale(
                "log" if event.inaxes.get_xscale() == "linear" else "linear"
            )
        elif key == "q":
            plt.close("all")
            sys.exit(0)

        if key == "w":
            self.ax.autoscale(enable=True, axis="both", tight=True)
            self.ax.relim()
            self.ax.autoscale_view()
            if hasattr(self, "ax_ratio") and self.ax_ratio.get_visible():
                self.ax_ratio.set_xlim(self.ax.get_xlim())
            self.draw_lines()
        elif key == "m":
            if event.inaxes == self.ax:
                self.identify_line(event.xdata)

        if self.continuum_mode and event.inaxes == self.ax:
            if key == "x":
                self.add_anchor(event.xdata, event.ydata)
            elif key == "d":
                self.delete_anchor(event.xdata, event.ydata)
            elif key == "c":
                self.update_fit()
            elif key == " ":
                self.handle_freeze(event.xdata)
            elif key == "u":
                self.handle_unfreeze(event.xdata)
            elif key == "a":
                self.handle_accept_all()

        if key == "e" and event.inaxes == self.ax:
            self.handle_ew(event.xdata)

        if key == "z":
            self.show_line_check_window()

        self.fig.canvas.draw_idle()

    def handle_freeze(self, x):
        if self.freeze_start is None:
            self.freeze_start = x
            self.freeze_guide = self.ax.axvline(
                x, color="grey", linestyle="--", alpha=0.5
            )
            print(f"Freeze start at {x:.2f}. Press Space again to end.")
        else:
            start = min(self.freeze_start, x)
            end = max(self.freeze_start, x)
            print(f"Freezing range: {start:.2f} - {end:.2f}")
            if self.freeze_guide:
                self.freeze_guide.remove()
                self.freeze_guide = None
            self.freeze_start = None

            if not self.continuum_fit_line or self.current_poly is None:
                print("No active fit to freeze.")
                return

            mask_indices = np.where(
                (self.wavelength >= start) & (self.wavelength <= end)
            )[0]
            if len(mask_indices) == 0:
                print("No data points in frozen range.")
                return

            y_values = self.current_poly(self.wavelength[mask_indices])
            self.frozen_segments.append(
                {"mask_idx": mask_indices, "y": y_values, "range": (start, end)}
            )
            self.update_continuum_display()

    def handle_unfreeze(self, x):
        to_remove = [
            i
            for i, seg in enumerate(self.frozen_segments)
            if seg["range"][0] <= x <= seg["range"][1]
        ]
        for i in sorted(to_remove, reverse=True):
            print(f"Unfreezing segment {i}")
            self.frozen_segments.pop(i)
        if to_remove:
            self.update_continuum_display()

    def handle_accept_all(self):
        if self.current_poly is None:
            print("No active fit to accept.")
            return

        coverage = np.zeros(len(self.wavelength), dtype=bool)
        for seg in self.frozen_segments:
            coverage[seg["mask_idx"]] = True

        gap_indices = np.where(~coverage)[0]
        if len(gap_indices) > 0:
            diffs = np.diff(gap_indices)
            breaks = np.where(diffs > 1)[0]
            starts = np.insert(breaks + 1, 0, 0)
            ends = np.append(breaks, len(gap_indices) - 1)
            added_count = 0
            for s, e in zip(starts, ends):
                indices = gap_indices[s : e + 1]
                if len(indices) == 0:
                    continue
                start_w = self.wavelength[indices[0]]
                end_w = self.wavelength[indices[-1]]
                y_vals = self.current_poly(self.wavelength[indices])
                self.frozen_segments.append(
                    {"mask_idx": indices, "y": y_vals, "range": (start_w, end_w)}
                )
                added_count += 1
            print(f"Accepted fit into {added_count} new segments.")
            self.update_continuum_display()

        self.model_accepted = True
        self.btn_save.label.set_text("Save Fit")
        self.btn_save.color = "0.85"
        self.btn_save.hovercolor = "0.95"
        print("Model Accepted. You can now save the fit.")
        self.fig.canvas.draw_idle()

    def update_continuum_display(self):
        if not self.frozen_segments:
            if self.frozen_line_obj:
                self.frozen_line_obj.remove()
                self.frozen_line_obj = None
            self.fig.canvas.draw_idle()
            return

        sum_y = np.zeros(len(self.wavelength))
        count_y = np.zeros(len(self.wavelength))
        for seg in self.frozen_segments:
            idx = seg["mask_idx"]
            sum_y[idx] += seg["y"]
            count_y[idx] += 1

        mask_valid = count_y > 0
        self.continuum_model = np.full_like(sum_y, np.nan)
        self.continuum_model[mask_valid] = sum_y[mask_valid] / count_y[mask_valid]

        if self.frozen_line_obj:
            self.frozen_line_obj.remove()

        (self.frozen_line_obj,) = self.ax.plot(
            self.wavelength,
            self.continuum_model,
            color="grey",
            linewidth=2,
            linestyle="-",
        )
        self.update_ratio_plot()

    def save_fit(self, event):
        if not self.model_accepted:
            print("Please press 'a' to accept/lock the full model first.")
            return

        if not hasattr(self, "continuum_model") or np.all(
            np.isnan(self.continuum_model)
        ):
            print("No valid continuum model found.")
            return

        base, ext = os.path.splitext(self.filename)
        cont_fname = f"{base}_cont{ext}" if ext else f"{base}_cont.txt"
        norm_fname = f"{base}_norm{ext}" if ext else f"{base}_norm.txt"

        is_anchor = np.zeros(len(self.wavelength), dtype=int)
        for ax, ay in self.continuum_anchors:
            idx = (np.abs(self.wavelength - ax)).argmin()
            is_anchor[idx] = 1

        try:
            np.savetxt(
                cont_fname,
                np.column_stack((self.wavelength, self.continuum_model, is_anchor)),
                header="IsContAnchor",
            )
            norm_flux = self.flux / self.continuum_model
            norm_err = self.error / self.continuum_model
            np.savetxt(
                norm_fname,
                np.column_stack((self.wavelength, norm_flux, norm_err, is_anchor)),
                header="Wavelength NormFlux NormError IsAnchor",
            )
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
  m         : Identify line at cursor
  z         : Open Line Check Window (3x2 grid)
  Redshift  : Enter value in textbox

Continuum Fitting:
  Continuum checkbox  : Enable mode + ratio panel
  Make Anchors        : Sigma-clipped median from flux (100px chunks)
  Fit Spline          : Fit spline through current anchors
  Auto Continuum      : Auto-detects <name>_cont<ext>, loads it,
                        places anchors on continuum values, and
                        opens the ratio panel automatically.
                        Button turns GREEN when the file exists.
  x  : Add anchor at cursor
  d  : Delete nearest anchor
  c  : Recompute spline fit

Frozen Segments:
  Space x2 : Freeze region (locks current fit there)
  u         : Unfreeze segment at cursor
  a         : Accept All → enables Save

Equivalent Width:
  e x2 : Define integration window
         Uses: loaded model > active poly > linear/1.0 fallback

Saving:
  Save Fit : Writes _cont and _norm files
             (enabled after Accept All or loading a continuum)
"""
        hfig = plt.figure(figsize=(6, 9))
        hfig.canvas.manager.set_window_title("Spectra Viewer Help")
        hfig.text(
            0.05,
            0.95,
            help_text,
            fontsize=10,
            family="monospace",
            verticalalignment="top",
        )
        plt.show()

    def handle_ew(self, x):
        if self.ew_start is None:
            self.ew_start = x
            self.ew_guide = self.ax.axvline(x, color="green", linestyle="--", alpha=0.5)
            print(f"EW Integration Start: {x:.2f}. Press 'e' again to end.")
            self.text_ew_obs.set_text("Obs: --")
            self.text_ew_rest.set_text("Rest: --")
            if self.ew_poly:
                self.ew_poly.remove()
                self.ew_poly = None
            self.fig.canvas.draw_idle()
        else:
            start = min(self.ew_start, x)
            end = max(self.ew_start, x)
            print(f"EW Integration Range: {start:.2f} - {end:.2f}")
            if self.ew_guide:
                self.ew_guide.remove()
                self.ew_guide = None
            self.ew_start = None

            mask = (self.wavelength >= start) & (self.wavelength <= end)
            if not np.any(mask):
                print("No data in range.")
                return

            wave_seg = self.wavelength[mask]
            flux_seg = self.flux[mask]

            if hasattr(self, "continuum_model") and not np.all(
                np.isnan(self.continuum_model[mask])
            ):
                cont_seg = self.continuum_model[mask]
                mode = "Global Model"
            elif self.current_poly is not None:
                cont_seg = self.current_poly(wave_seg)
                mode = "Active Poly"
            else:
                if self.linear_cont_fallback:
                    idx_s = np.abs(self.wavelength - start).argmin()
                    idx_e = np.abs(self.wavelength - end).argmin()
                    y1, x1 = self.flux[idx_s], self.wavelength[idx_s]
                    y2, x2 = self.flux[idx_e], self.wavelength[idx_e]
                    gradient = (y2 - y1) / (x2 - x1) if x2 != x1 else 0
                    cont_seg = y1 + gradient * (wave_seg - x1)
                    mode = "Linear Fallback"
                else:
                    cont_seg = np.ones_like(flux_seg)
                    mode = "Normalized (1.0)"

            ew_obs = np.trapz(1.0 - flux_seg / cont_seg, x=wave_seg)
            self.ew_obs_value = ew_obs
            ew_rest = ew_obs / (1 + self.current_z)

            print(f"EW ({mode}): Obs={ew_obs:.3f}, Rest={ew_rest:.3f}")
            self.text_ew_obs.set_text(f"Obs: {ew_obs:.3f} A")
            self.text_ew_rest.set_text(f"Rest: {ew_rest:.3f} A")

            self.ew_poly = self.ax.fill_between(
                wave_seg, cont_seg, flux_seg, color="green", alpha=0.3
            )
            self.fig.canvas.draw_idle()

    def auto_continuum(self, event):
        if not self.continuum_mode:
            print("Please enable Continuum Mode first.")
            return

        for pt in self.continuum_plots:
            pt.remove()
        self.continuum_anchors = []
        self.continuum_plots = []
        print("Auto-generating anchors from flux...")

        chunk_size = 100
        n_pixels = len(self.wavelength)

        for i in range(0, n_pixels, chunk_size):
            end = min(i + chunk_size, n_pixels)
            if end - i < 10:
                continue

            wave_chunk = self.wavelength[i:end]
            flux_chunk = self.flux[i:end]
            mean, median, std = sigma_clipped_stats(flux_chunk, sigma=3.0, maxiters=5)

            if np.ma.is_masked(median):
                continue

            anchor_x = float(np.mean(wave_chunk))
            anchor_y = float(median)

            in_frozen = any(
                seg["range"][0] <= anchor_x <= seg["range"][1]
                for seg in self.frozen_segments
            )
            if not in_frozen:
                self.add_anchor(anchor_x, anchor_y)

        print(f"Added {len(self.continuum_anchors)} anchors.")

    def add_anchor(self, x, y):
        for seg in self.frozen_segments:
            if seg["range"][0] <= x <= seg["range"][1]:
                print("Cannot add anchor in frozen region!")
                return
        self.continuum_anchors.append((x, y))
        (pt,) = self.ax.plot(x, y, "m*", markersize=10)
        self.continuum_plots.append(pt)
        self.fig.canvas.draw_idle()

    def delete_anchor(self, x, y):
        if not self.continuum_anchors:
            return
        dists = [
            np.sqrt((ax - x) ** 2 + (ay - y) ** 2) for ax, ay in self.continuum_anchors
        ]
        idx = np.argmin(dists)
        anchor = self.continuum_anchors[idx]
        for seg in self.frozen_segments:
            if seg["range"][0] <= anchor[0] <= seg["range"][1]:
                print("Cannot delete anchor in frozen region!")
                return
        self.continuum_anchors.pop(idx)
        pt = self.continuum_plots.pop(idx)
        pt.remove()
        self.fig.canvas.draw_idle()

    def update_fit(self):
        active_anchors = [
            p
            for p in self.continuum_anchors
            if not any(
                seg["range"][0] <= p[0] <= seg["range"][1]
                for seg in self.frozen_segments
            )
        ]
        if len(active_anchors) < self.poly_order + 1:
            print(f"Not enough active anchors for order {self.poly_order} fit.")
            return

        anchors = sorted(active_anchors, key=lambda p: p[0])
        x_anchors = np.array([p[0] for p in anchors])
        y_anchors = np.array([p[1] for p in anchors])

        try:
            self.current_poly = UnivariateSpline(
                x_anchors, y_anchors, k=self.poly_order, s=0
            )
            xmin, xmax = self.ax.get_xlim()
            x_plot = np.linspace(xmin, xmax, 1000)
            y_plot = self.current_poly(x_plot)

            if self.continuum_fit_line:
                self.continuum_fit_line.remove()

            (self.continuum_fit_line,) = self.ax.plot(
                x_plot,
                y_plot,
                color="magenta",
                linestyle="dashed",
                linewidth=2,
                label="Spline Fit",
            )
            self.update_ratio_plot()
            print("Fitted spline to anchors.")
        except Exception as e:
            print(f"Fit failed: {e}")

    def fit_spline(self, event):
        active_anchors = [
            p
            for p in self.continuum_anchors
            if not any(
                seg["range"][0] <= p[0] <= seg["range"][1]
                for seg in self.frozen_segments
            )
        ]
        if len(active_anchors) < 2:
            print("Not enough active anchors for spline fit (need at least 2).")
            return

        anchors = sorted(active_anchors, key=lambda p: p[0])
        x_anchors = np.array([p[0] for p in anchors])
        y_anchors = np.array([p[1] for p in anchors])

        try:
            self.current_poly = UnivariateSpline(
                x_anchors, y_anchors, k=self.poly_order, s=0
            )
            xmin, xmax = self.ax.get_xlim()
            x_plot = np.linspace(xmin, xmax, 1000)
            y_plot = self.current_poly(x_plot)

            if self.continuum_fit_line:
                self.continuum_fit_line.remove()

            (self.continuum_fit_line,) = self.ax.plot(
                x_plot,
                y_plot,
                color="magenta",
                linestyle="dashed",
                linewidth=2,
                label="Spline Fit",
            )
            self.update_ratio_plot()
            print("Fitted spline to anchors.")
        except Exception as e:
            print(f"Spline fit failed: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python spectra_viewer.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    viewer = SpectraViewer(filename)
