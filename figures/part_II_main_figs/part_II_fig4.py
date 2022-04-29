from matplotlib import pyplot as plt
from pyOER import StabilityResultsCollection


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


if False:  # set this to False after first time to run faster by loading collection.json
    # The names of each type of sample.
    sample_mapping = {
        "RT-RuO2": ["Taiwan", "Easter"],
        "400C-RuO2": ["Maundy", "Stoff", "Sofie", "Mette", "John"],
        "cycled RuO2": [
            "Taiwan1G",
        ],
        "RuOx/Ru": [
            "Bernie",
        ],
        "Ru foam": [
            "Evans12",
        ],
        "RT-IrO2": ["Goof", "Decade"],
        "400C-IrO2": ["Legend", "Champ"],
        "cycled IrO2": ["Decade1G", "Legend4C"],
        "IrOx/Ir": [
            "Jazz",
        ],
    }
    current_point_mapping = {
        s: float(s) for s in ["1.0", "0.5", "0.25", "0.15", "0.05"]
    }

    results_collection = StabilityResultsCollection(
        sample_mapping=sample_mapping, current_point_mapping=current_point_mapping
    )
    results_collection.save("results_collection.json")
else:  # load tof collection
    results_collection = StabilityResultsCollection.open("results_collection.json")

sample_plot_specs = {
    "RT-RuO2": {"color": "#54bdebff", "marker": "o"},
    "400C-RuO2": {"color": "#163854ff", "marker": "o"},
    "cycled RuO2": {"color": "b", "marker": "*", "markersize": 8},
    "RuOx/Ru": {"color": "b", "marker": "s"},
    "Ru foam": {"color": "m", "marker": "s"},
    "RT-IrO2": {"color": "#54ebbdff", "marker": "o"},
    "400C-IrO2": {"color": "#165438ff", "marker": "o"},
    "cycled IrO2": {"color": "g", "marker": "*", "markersize": 8},
    "IrOx/Ir": {"color": "g", "marker": "s"},
}

fig3b, ax3b = plt.subplots()
stats_collection = {}
fig3b_current_point = "0.5"
for sample_type in results_collection:
    specs = sample_plot_specs[sample_type]
    for tof_time in ["start", "steady"]:

        specs["markerfacecolor"] = "w" if tof_time == "start" else specs["color"]

        stats = results_collection.get_stats(
            sample_type=sample_type,
            current_point=fig3b_current_point,
            tof_time=tof_time,
        )

        S_number, sigma_S = stats["S_number"]
        S_number_lattice, sigma_S_lattice = stats["S_number_lattice"]

        ax3b.plot(S_number_lattice, S_number, **specs)
        if sigma_S:
            diss_error_specs = specs.copy()
            diss_error_specs.update(marker="_", linestyle="-", markersize=5)
            ax3b.plot(
                [S_number_lattice, S_number_lattice],
                [S_number - sigma_S, S_number + sigma_S],
                **diss_error_specs,
            )
        if sigma_S_lattice:
            exc_error_specs = specs.copy()
            exc_error_specs.update(marker="|", linestyle="-", markersize=5)
            ax3b.plot(
                [
                    S_number_lattice - sigma_S_lattice,
                    S_number_lattice + sigma_S_lattice,
                ],
                [S_number, S_number],
                **exc_error_specs,
            )

lims = [5, 1e5]
ax3b.set_xlim(lims)
ax3b.set_ylim(lims)
ax3b.plot(lims, lims, "k--")

ax3b.set_xlabel("oxygen stability number")
ax3b.set_ylabel("metal stability number")
ax3b.set_xscale("log")
ax3b.set_yscale("log")

if forpublication:
    fig3b.savefig("paper_II_v4_fig4.png")
    fig3b.savefig("paper_II_v4_fig4.svg")
