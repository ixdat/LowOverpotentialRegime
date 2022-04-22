from matplotlib import pyplot as plt

from pyOER import Experiment


forpublication = True
if forpublication:  # for the publication figure
    import matplotlib as mpl

    mpl.rcParams["figure.figsize"] = (3.25, 2.75)
    # plt.rc('text', usetex=True)  # crashingly slow
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=8)
    plt.rc("lines", linewidth=0.6)
else:
    plt.style.use("default")


if True:  # fig 2a, Reshma1 in 16-O electrolyte
    exp = Experiment.open(47)

    if True:  # Faradaic efficiency plot
        axes = exp.plot_faradaic_efficiency()
        fig = axes[0].get_figure()
        fig.subplots_adjust(right=0.85)
        if not forpublication:
            axes[1].set_title(str(exp))

    if True:  # activity plot (fig 2a)
        # exp.correct_current()
        exp.measurement.meas.reset_bg()  # so that background is not subtracted
        exp.measurement.cut_meas(tspan=[1700, 5400], t_zero="start")
        axes = exp.measurement.plot(
            mol_list=list(exp.mdict.values()),
            unit="pmol/s/cm^2",
            removebackground=False,
        )
        axes[0].set_ylabel("O$_2$ / (pmol s$^{-1}$cm$^{-2}_{geo})$")
        axes[1].set_ylabel("E vs RHE / (V)")
        axes[3].set_ylabel("J / (mA cm$^{-2}_{geo}$)")
        axes[0].set_xlabel("time / (s)")
        axes[1].set_xlabel("time / (s)")
        axes[0].set_yscale("log")
        axes[0].set_ylim([5e-1, 3e3])
        if forpublication:
            axes[0].get_figure().savefig(f"paper_I_v5_fig2a.png")
            axes[0].get_figure().savefig(f"paper_I_v5_fig2a.svg")
        else:
            axes[1].set_title(str(exp))


if True:  # fig 2b, Reshma1 in 18-O electrolyte  (wierd FE)
    exp = Experiment.open(55)

    if True:  # Faradaic efficiency plot
        axes = exp.plot_faradaic_efficiency()
        fig = axes[0].get_figure()
        fig.subplots_adjust(right=0.85)
        if forpublication:
            axes[1].set_title(str(exp))

    if True:  # activity plot (fig 2a)
        # exp.correct_current()
        exp.measurement.meas.reset_bg()  # so that background is not subtracted
        exp.measurement.cut_meas(tspan=exp.tspan_plot, t_zero="start")
        axes = exp.measurement.plot(
            mol_list=list(exp.mdict.values()),
            unit="pmol/s/cm^2",
            removebackground=False,
        )
        axes[0].set_ylabel("O$_2$ / (pmol s$^{-1}$cm$^{-2}_{geo})$")
        axes[1].set_ylabel("E vs RHE / (V)")
        axes[3].set_ylabel("J / (mA cm$^{-2}_{geo}$)")
        axes[0].set_xlabel("time / (s)")
        axes[1].set_xlabel("time / (s)")
        axes[0].set_yscale("log")
        axes[0].set_ylim([5e-1, 3e3])
        if forpublication:
            axes[0].get_figure().savefig(f"paper_I_v5_fig2b.png")
            axes[0].get_figure().savefig(f"paper_I_v5_fig2b.svg")
        else:
            axes[1].set_title(str(exp))

if True:  # fig 2c, Reshma1 in 18-O electrolyte  (sensible FE)
    exp = Experiment.open(54)
    exp.correct_current()

    if True:  # Faradaic efficiency plot
        axes = exp.plot_faradaic_efficiency()
        fig = axes[0].get_figure()
        fig.subplots_adjust(right=0.85)
        axes[1].set_ylim([0, 0.32])
        axes[0].set_ylabel("Faradaic efficiency / (%)")
        axes[0].set_xlabel("E vs RHE / (V)")
        axes[1].set_ylabel("Current density / (mA cm$^{-2}_{geo}$)")
        if forpublication:
            fig.savefig(f"paper_I_v5_fig2c.png")
            fig.savefig(f"paper_I_v5_fig2c.svg")
        else:
            axes[1].set_title(str(exp))

    if True:  # activity plot (fig 2a)
        exp.measurement.meas.reset_bg()
        exp.measurement.cut_meas(tspan=exp.tspan_plot, t_zero="start")
        axes = exp.measurement.plot(
            mol_list=list(exp.mdict.values()),
            unit="pmol/s/cm^2",
            removebackground=False,
        )
        if not forpublication:
            axes[1].set_title(str(exp))

if True:  # fig 2d, with the inset
    exp = Experiment.open(71)
    exp.measurement.plot(
        mol_list=[exp.mdict["O2_M36"], exp.mdict["O2_M34"]],
        unit="pmol/s/cm^2",
        removebackground=True,
        tspan=[0, 5000],
        logplot=False,
    )  # garbage plot
    exp.measurement.meas.reset_bg()
    axes = exp.measurement.plot(
        mol_list=[exp.mdict["O2_M34"], exp.mdict["O2_M36"]],
        unit="pmol/s/cm^2",
        removebackground=False,
        tspan=[0, 5000],
        logplot=False,
    )  # garbage plot
    axes[0].set_ylabel("O$_2$ / (pmol s$^{-1}$cm$^{-2}_{geo})$")
    axes[1].set_ylabel("E vs RHE / (V)")
    axes[3].set_ylabel("J / (mA cm$^{-2}_{geo}$)")
    axes[0].set_xlabel("time / (s)")
    axes[1].set_xlabel("time / (s)")
    axes[0].set_yscale("log")
    axes[0].set_ylim([5, 3e3])
    if forpublication:
        axes[0].get_figure().savefig(f"paper_I_v5_fig2d.png")
        axes[0].get_figure().savefig(f"paper_I_v5_fig2d.svg")
    else:
        axes[1].set_title(str(exp))

    if True:  # inset
        tspan_inset = [2600, 4000]
        x, y = exp.measurement.meas.grab_flux(
            exp.mdict["O2_M36"],
            tspan=tspan_inset,
            removebackground=False,
        )
        y *= 1e12 / exp.meas.A_el  # [mol/s] -> [pmol/s/cm^2]
        t, V = exp.measurement.meas.grab("potential", tspan=tspan_inset)

        fig, ax = plt.subplots()
        ax2 = ax.twinx()

        if False:  # fade out original data, include smoothed data

            from EC_MS import smooth

            n_points = 15
            y_smooth = smooth(y, n_points=n_points)  # 10-point moving average
            ax.plot(x, y, "g", alpha=0.2)
            ax.plot(x[n_points:-n_points], y_smooth[n_points:-n_points], "g")
        else:  # plot original data without fading.
            ax.plot(x, y, "g")
            ax.set_ylim([2, 7])
        # ax.set_ylim([3.3, 5.8])
        ax2.plot(t, V, "k")
        ax2.set_ylim(1.19, 1.65)
        ax.set_yticks([])
        ax2.set_xticks([])
        ax2.set_yticks([])

        if forpublication:
            fig.set_figwidth(1)
            fig.set_figheight(1)
            fig.savefig(f"paper_I_v6_fig2d_inset.png")
            fig.savefig(f"paper_I_v6_fig2d_inset.svg")
