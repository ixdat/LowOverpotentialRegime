import pyOER
import datetime
import numpy as np
import matplotlib.pyplot as plt

forpublication = True
tex = False
if forpublication:  # for the publication figure
    import matplotlib as mpl

    mpl.rcParams["figure.figsize"] = (3.25, 2.75)
    if tex:
        plt.rc('text', usetex=tex)  # crashingly slow
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=8)
    plt.rc("lines", linewidth=0.6)
    plt.rc("lines", markersize=3) # extra
else:
    plt.style.use("default")

sputter_estimation = {
    'sample': 'Easter1A',
    'keys': {
        # Start sputtering in a new spot on the sample (before measurement)
        # active: label
        7: 'New spot',
        15: '5min He',
        12: '30s Ar',
        13: '60s Ar',
        11: '90s Ar',
        14: '120s Ar',
        9: '180s Ar',
        10: '300s Ar',
        6: '600s Ar',
    },
}

dict_of_samples = [{
    'Decade1A': {# m231
        'before': 2,
        'after': 0,
        'sputter': 1,
        'type': 'RT-IrO2',
        'isotope': 18,
        'notes': '',
        }},{
    'Decade1B': {# m235
        'before': 1,
        'after': None,
        'sputter': None,
        'type': 'RT-IrO2',
        'isotope': 18,
        'notes': '',
        }},{
    'Decade1G': {# m233
        'before': 7,
        'after': 1,
        'sputter': 3,
        'type': 'cycled IrO2',
        'isotope': 18,
        'notes': '',
        }},{

    # Used for sputter rate estimation
    'Easter1A': {# m160
        'before': 5,
        'after': 4, # maybe 0 !!
        # Unclear from data collection
        'sputter': 0, # 2: 30s and 3: 60s and 1: "point-2"
        'type': 'RT-RuO2',
        'isotope': 18,
        'notes': 'tested for sputter rate. Measurements AFTER ambiguous.',
        }},{

    # Also used for sputter rates
    'Easter1B': {# m197, m198
        'before': 14,
        'after': None, # 0
        'sputter': None, # 1
        'type': 'RT-RuO2',
        'isotope': 18,
        'notes': 'm[168,170,161] before LEIS. Heavily corroded: AFTER/SPUTTER disregarded',
        }},{

    'Easter1C': {# m195
        'before': 4, # 0 is Goof1 (A)
        'after': None, # 1
        'sputter': None, # 2
        'type': 'RT-RuO2',
        'isotope': 18,
        'notes': 'Active 0 is really Goof1(A) BEFORE. Heavily corroded: AFTER/SPUTTER disregarded',
        }},{

    'Folk3': {# m164, m203
        'before': 0,
        'after': 3,
        'sputter': 2,
        'type': 'EC 400C-IrO2', # TODO correct?
        'isotope': 16,
        'notes': '', # Check notes: Goof1 / Taiwan_after_first_testing
        }},{

    'Goof1A': {# m165, m169
        'before': None, # 0 in Easter1C
        'after': 0,
        'sputter': None,
        'type': 'RT-IrO2',
        'isotope': 18,
        'notes': '',
        }},{
        
    'Goof1B': {# m171
        'before': 3,
        'after': None,
        'sputter': None,
        'type': 'RT-IrO2',
        'isotope': 18,
        'notes': '2 is void. 1 and 0 are clearly labelled as 16-O (10 days later)',
        },
        
    'John4A': {# m222
        'before': 2, # check 3
        'after': 1,
        'sputter': 0,
        'type': '400C-RuO2',
        'isotope': 18,
        'notes': '', # Note: 19L22_John_and_Love
        }},{
        
    'John4C': {# m223
        'before': 1,
        'after': 0,
        'sputter': None,
        'type': '400C-RuO2',
        'isotope': 18,
        'notes': '', # Note: 19L22_John_and_Love
        }},{
        
    'Legend4A': {# m224
        'before': 0,
        'after': None,
        'sputter': None,
        'type': '400C-IrO2',
        'isotope': 18,
        'notes': '',
        }},{

    'Legend4C': {# m227
        'before': 0,
        'after': 6, # 5 (corroded)
        'sputter': 8, # 7 (corroded)
        'type': 'cycled IrO2',
        'isotope': 18,
        'notes': 'outlier'# What is m227?
        }},{

    'Legend4C': {# m232
        'before': 8, # presputtered
        'after': 2, # 1 (corroded)
        'sputter': 4, # 3 (corroded)
        'type': 'cycled IrO2',
        'isotope': 18,
        'notes': 'Used for m227 before m232 (not a clean experiment).',
        }},{

    'Legend4C': {# m232
        'before': 8, # presputtered
        'after': 1, # 1 (corroded)
        'sputter': 3, # 3 (corroded)
        'type': 'cycled IrO2',
        'isotope': 18,
        'notes': 'Used for m227 before m232 (not a clean experiment).',
        }},{

    # Special sample
    'Melih2': {# m157 m158
        'before': None,#(3, 4),
        'after': None,#0,
        'sputter': None,
        'type': 'UHV 18O',
        'isotope': 18,
        'notes': 'Leis AFTER is long after measurement! Also 5d contaminant.',
        }},{

    'Nancy1': {# m84 m99
        'before': None,
        'after': (0, 2), # (0, 6) for high res, but sputtered more with He
        'sputter': (1, 2),
        'type': 'EC RT-RuO2', # TODO correct?
        'isotope': 16,
        'notes': '',
        }},{

    'Nancy5': {# m89 m102
        'before': None,
        'after': (1, 2),
        'sputter': (0, 1),
        'type': 'EC RT-RuO2', # TODO correct?
        'isotope': 16,
        'notes': '',
        }},{

    'Nancy6': {# m90 m91 m92 m103
        'before': None,
        'after': (1, 2),
        'sputter': 0,
        'type': 'EC RT-RuO2', # TODO correct?
        'isotope': 16,
        'notes': 'Titanium present in small amount.',
        }},{

    'Reshma4E': {# m159
        'before': 1,
        'after': 2,
        'sputter': None,
        'type': '400C-RuO2',
        'isotope': 16,
        'notes': '',
        }},{

    'Reshma4F': {# m187
        'before': 2,
        'after': 0,
        'sputter': 1,
        'type': '400C-RuO2',
        'isotope': 16,
        'notes': '',
        }},{

    'Reshma4I': {# m234
        'before': None,
        'after': 0,
        'sputter': 1,
        'type': '400C-RuO2',
        'isotope': 16,
        'notes': '',
        }},{

    'Stoff4A': {# m192
        'before': 4,
        'after': 0,
        'sputter': 1, # much higher signal than the other spectra
        'type': '400C-RuO2',
        'isotope': 18,
        'notes': '',
        }},{

    'Stoff4B': {# m193
        'before': 2, # also 3 at lower signal. The isotope ratio may have changed sign.
        'after': 0,
        'sputter': 1,
        'type': '400C-RuO2',
        'isotope': 18,
        'notes': '',
        }},{

    'Taiwan1A': {# m204
        'before': 1,
        'after': 3,
        'sputter': 2,
        'type': 'RT-RuO2',
        'isotope': 18,
        'notes': '', # What is active=4?
        }},{

    'Taiwan1C': {# m207 m208 m209
        'before': 1,
        'after': 4,
        'sputter': 2,
        'type': 'RT-RuO2',
        'isotope': 18,
        'notes': '',
        }},{

    'Taiwan1D': {# m230
        'before': 4, # Change to None maybe
        'after': 0,
        'sputter': 1,
        'type': 'RT-RuO2',
        'isotope': 18,
        'notes': 'Long time between BEFORE and measurement',
        }},{

#    # Dark region
#    'Taiwan1G': {# m218 m216 m217
#        'before': 6,
#        'after': 3,
#        'sputter': None,#0, # Lots of Ti resembling spectrum 2
#        'type': 'cycled RuO2',
#        'isotope': 18,
#        'notes': 'Dark region. SPUTTER disregarded because of titanium. Measured 3 different pos. Possibly mislabeled.',
#        }},{
#
#    # Light region
#    'Taiwan1G': {# m218 m216 m217
#        'before': 6,
#        'after': None,#4, # Lots of Ti
#        'sputter': 1, # Not even a trace of Ti
#        'type': 'cycled RuO2',
#        'isotope': 18,
#        'notes': 'Light region. SPUTTER disregarded because of titanium. Measured 3 different pos. Possibly mislabeled.',
#        }},{
#
#    # Center region
#    'Taiwan1G': {# m218 m216 m217
#        'before': 6,
#        'after': None,#5,
#        'sputter': None#2, # None (corrosion) Lots of Ti
#        'type': 'cycled RuO2',
#        'isotope': 18,
#        'notes': 'Center region. SPUTTER disregarded because of titanium. Measured 3 different pos. Possibly mislabeled.',
#        }},
    # Custom (combined) region
    'Taiwan1G': {# m218 m216 m217
        'before': 6,
        'after': 3,
        'sputter': 1,#0, # Lots of Ti resembling spectrum 2
        'type': 'cycled RuO2',
        'isotope': 18,
        'notes': 'Custom region. SPUTTER disregarded because of titanium. Measured 3 different pos. Possibly mislabeled.',
        }},
]

sample_plot_specs = {
    # Extra
    'EC RT-RuO2':   {"color": "c",          "marker": "d"},
    'EC 400C-IrO2': {"color": "y",          "marker": "d"},
    'UHV 18O':      {"color": "k",          "marker": "d"},# 'markerfacecolor': 'w'},

    "400C-RuO2":    {"color": "#163854ff",  "marker": "o"},
    "cycled RuO2":  {"color": "b",          "marker": "*", "markersize": 8},
    "RuOx/Ru":      {"color": "b",          "marker": "s"},
    "Ru foam":      {"color": "m",          "marker": "s"},
    "RT-IrO2":      {"color": "#54ebbdff",  "marker": "o"},
    "400C-IrO2":    {"color": "#165438ff",  "marker": "o"},
    "cycled IrO2":  {"color": "g",          "marker": "*", "markersize": 8},
    "IrOx/Ir":      {"color": "g",          "marker": "s"},
    "RT-RuO2":      {"color": "#54bdebff",  "marker": "o"},
}

# Load ISS handler
handle = pyOER.ISS()
interactive = False

if False:
    # This section used to manually go through and check data for each sample
    sample = "Folk3"
    handle.get_sample(sample)
    handle.verbose = True
    before = []
    after = []
    # Plot "Before" and "After"
    handle.plot(before, mass_lines=[16, 18, 35, 47.8])
    print('_'*20)
    handle.plot(after, mass_lines=[16, 18, 35, 47.8])

fig_individual = plt.figure('Individual results')
ax_individual = plt.gca()#fig_individual.add_axes(111)

fig_total = plt.figure('Total results')
ax_total = plt.gca()

xticks = [0, 1, 2]#-1, -2, -3]
xlabels = [
    'Before OER',
    'After OER',
    'OER \& Sputtered',
    ] if tex else [
    'Before OER',
    'After OER',
    'OER & Sputtered',
    ]

for ax in [ax_individual, ax_total]:
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=15)
    ax.set_yticks([0, 25, 50, 75, 100])
    ylabel = r'$^{16}$O content / ($\%$)' if tex else '$^{16}$O content / (%)'
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 80)
    ax.set_xlim(-0.1, 2.4)
    #ax_individual.set_ylim(0, 100)
if interactive:
    plt.ion()
    plt.show()
#1/0
total_results = {
    stype: {
        isotope: {
            'before': [],
            'after': [],
            'sputter': [],
            'x_mod': np.random.random() / 3,
            } for isotope in [16, 18]
        }
    for stype
    in sample_plot_specs
    }

for sample_dict in dict_of_samples:
    sample, info = next(iter(sample_dict.items()))
    if sample != handle.sample:
        handle.get_sample(sample)
    print(f' *** {sample} ***')
    specs = sample_plot_specs[info['type']]

    # Plot results (individual fits)
    if False:
        if 'isotope' in info:
            continue
        for order in ['before', 'after', 'sputter']:
            index = info[order]
            if index is not None:
                if isinstance(index, tuple):
                    handle.active = index[0]
                    handle.data.default_scan = index[1]
                elif isinstance(index, int):
                    handle.active = index
                else:
                    raise TypeError(f'Index "{index}" not understood')

                print(f'Plotting fit for {sample}: "{order}"')
                handle.plot_fit()
            else:
                print(f'No data of type "{order}"')

    if sample in ['Melih2']:
        continue
    if info['isotope'] == 16:
        continue
    if info['notes'] == 'outlier':
        continue
    # Plot results (all samples unto results graph)
    O16 = np.array([np.nan]*3)
    #x_mod = np.random.random() / 3
    for i, order in enumerate(['before', 'after', 'sputter']):
        index = info[order]
        if index is not None:
            if isinstance(index, tuple):
                handle.active = index[0]
                handle.data.default_scan = index[1]
            elif isinstance(index, int):
                handle.active = index
            else:
                raise TypeError(f'Index "{index}" not understood')

            print(f'Plotting fit for {sample}: "{order}"')
            O16[i] = handle.meta('O16')
            total_results[info['type']][info['isotope']][order].append(O16[i])
        else:
            print(f'No data of type "{order}"')
        x_mod = total_results[info['type']][info['isotope']]['x_mod']
        mask = np.isfinite(O16)
        x = (np.array(xticks) + x_mod)[mask]
        y = O16[mask]
        ax_individual.plot(
            x,
            y,
            linestyle='solid' if len(x) == 3 else 'dashed',
            mfc='white' if info['isotope'] == 16 else None,
            label=sample,
            **specs,
            )
        if interactive:
            ax_individual.set_title(
                f'{sample} {handle.active, handle.data.default_scan} - {O16[i]}',
                color=specs['color'],
                )
            if True:#sample == 'Folk3':
                print(f'\n{sample}: {handle.active, handle.data.default_scan}\n{order}')
                print('O16: ', O16)
                print('x: ', x)
                print('y: ', y)
                req = input('Plot fit? ')
                if req == 'p':
                    handle.plot_fit()
                    input('<Enter>')

for stype, specs in sample_plot_specs.items():
    #x_mod = np.random.random() / 3
    for isotope in [16, 18]:
        mean, std = [], []
        x_mod = total_results[stype][isotope]['x_mod']
        full_array = [int(bool(total_results[stype][isotope][order])) for order in ['before', 'after', 'sputter']]
        print(full_array)
        if sum(full_array) < 2:
            continue
        for i, order in enumerate(['before', 'after', 'sputter']):
            array = total_results[stype][isotope][order]
            print(array)
            if array:
                array = np.array(array)
                mean.append(np.average(array))
                std.append(np.std(array))
            else:
                mean.append(np.nan)
                std.append(np.nan)
            ax_total.errorbar(
                xticks[i] + x_mod,
                mean[i],
                yerr=std[i],
                markerfacecolor= 'white' if len(array) == 1 else specs['color'],
                **specs,
                )
            """
            ax_total.plot(
                xticks[i] + x_mod,
                mean[i],
                'o',
                ms=10,
                mew=2,
                mfc=d,
                mec=specs['color'],
                #**specs,
                )
            """
        mask = np.isfinite(mean)
        x = (np.array(xticks) + x_mod)[mask]
        y = np.array(mean)[mask]
        ax_total.plot(
            x,
            y,
            linestyle='solid' if len(x) == 3 else 'dashed',
            color=specs['color'],
            )

for ax in [ax_individual, ax_total]:
    plt.sca(ax)
    plt.tight_layout()

if True:
    plt.sca(ax_individual)
    name = 'individual_leis_results'
    plt.savefig(name + '.png')
    plt.savefig(name + '.svg')

    plt.sca(ax_total)
    name = 'total_leis_results'
    plt.savefig(name + '.png')
    plt.savefig(name + '.svg')

if interactive:
    import time
    time.sleep(5)
    input('End...')

plt.show()
