import numpy as np
import matplotlib.pyplot as plt
import pyOER

samples = {}
# RuO2
samples["RuO2 amorphous"] = [
    "Nancy1",  # 16
    "Nancy5",  # 16
    "Nancy6",  # 16
    "Easter1A",  #    18
    "Easter1B",  #    18
    "Easter1C",  #    18
    "Taiwan1A",  #    18
    "Taiwan1C",  #    18
    "Taiwan1G",  #    18
    "Taiwan1D",  #    18
]
samples["RuO2 rutile"] = [
    "Reshma4E",  # 16 (400 C)
    "Reshma4F",  # 16 (400 C)
    "Reshma4I",  # 16 (400 C)
    "Stoff4E",  #    18
    "Stoff4F",  #    18
    "Stoff4A",  #    18
    "Stoff4B",  #    18
    "Stoff1D",  #    18
    "John4A",  #    18
    "John4C",  #    18
]
samples["Ru foam"] = [
    "Evans9",  # 16
    "Evans7",  # 16
    "Evans2",  # 16
    "Evans8",  # 16
    "Evans12",  # 16
]
samples["Ru metallic"] = [
    "Bernie4",  # EC-treated
    "Bernie5",  # EC-treated
    "Melih2",  # treated
]
# Pt
samples["Pt"] = [
    "Trimi1",
]
# Ir (all)
samples["Ir"] = [
    "Jazz5",  # EC-treated (metallic)
    "Folk3",  # EC-treated (rutile)
    "Goof1A",  # RT 18
    "Goof1B",  # RT 18
    "Legend4A",  # 400C 18
    "Legend4C",  # 400C 18
    "Decade1A",  # RT 18
    "Decade1G",  # RT 18
    "Decade1B",  # RT 18
    "Decade1C",  # RT 18
]

# Choose selected groups here ### SELECTION ###
selection = [
    #'RuO2 amorphous',
    "RuO2 rutile",
    #'Ru foam',
    #'Ru metallic',
    #'Pt',
    #'Ir',
]
title = " _ ".join(selection)

colors = ["k", "r", "g", "b", "m"] * 10
invalid_samples = []
datas = {}
ratios = {}
names = [sample_ for selected in selection for sample_ in samples[selected]]

# Calculate ratios
for sample_counter, sample in enumerate(names):
    try:
        data = pyOER.ISS(sample)
    except FileNotFoundError:
        print(f"Sample {sample} not found...")
        continue
    if len(data.keys) == 0:
        print(f'Could not find data matching "{sample}"')
        invalid_samples.append(sample_counter)
        continue
    datas[sample] = data
    print("*" * 20)
    print(f"Available keys: {data.keys}")
    ratios[sample], coeffs = data.fit_with_reference(
        peaks=[[16, 18]],
        plot_result=False,
    )
    for i in ratios[sample].keys():
        if data._active[i].good is False:
            continue
        print(data._active[i].filename)
        print(data._active[i].date)
        print(f'\nOxygen 16 content: {data.fit_ratios[i]["16"]*100} %\n')
        if False:
            # Visualize the resulting fits
            data.plot_fit(i, show=False)

# Remove invalid samples from list and inform user
invalid_samples.reverse()
for i in invalid_samples:
    sample = names.pop(i)
    print(f'Could not find data matching "{sample}"')

# Plot all O-16 ratios
plot_data = []
counter = 0
fig = plt.figure(title)
ax = fig.add_axes([0.05, 0.15, 0.9, 0.6])
for j, sample in enumerate(names):
    for i in datas[sample].keys:
        # Skip bad data
        if datas[sample]._active[i].good is False:
            continue
        # Plot good data
        plt.plot(counter, ratios[sample][i]["16"] * 100, "o", color=colors[j])
        plot_data.append(
            [
                sample,
                datas[sample],
                datas[sample]._active[i].sample,
                datas[sample]._active[i].date,
                counter,
                ratios[sample][i]["16"],
                ratios[sample][i]["18"],
            ]
        )
        counter += 1
xticks = [i for (gen_name, data_object, name, date, i, r1, r2) in plot_data]
dates = [
    pyOER.iss.date_formatter(date)
    for (gen_name, data_object, name, date, i, r1, r2) in plot_data
]
xlabels = [
    f"{gen_name} {name.lstrip(gen_name)}"
    for (gen_name, data_object, name, date, i, r1, r2) in plot_data
]

secaxx = ax.secondary_xaxis("top")
secaxy = ax.secondary_yaxis("right")

# Update canvas
fig.canvas.draw()

secaxy.set_ylabel("O-18 ratio (%)")
yticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
ax.set_yticks(yticks)
ax.set_yticklabels(yticks)
secaxy.set_yticks(yticks)
yticks.reverse()
secaxy.set_yticklabels(yticks)
secaxx.set_xticks(xticks)
secaxx.set_xticklabels(dates, rotation=90, fontsize=12)
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, rotation=90, fontsize=12)
ax.set_ylabel("O-16 ratio (%)")
plt.grid(True)
plt.show()
