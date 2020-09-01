# @Author: Tristan Croll <tic20>
# @Date:   29-Nov-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 15-May-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



from math import radians, pi
import os
import matplotlib
matplotlib.rc('font', **{'family':'serif','serif':['FreeSerif']})

# When mutating, some torsions require adjustment to achieve the fairest overlap
# between before and after. In threonine, for example, chi1 is defined as
# N-CA-CB-OG1, which puts it 120 degrees out from valine. In general, it's best
# to rotate to give the longest possible overlap of carbon atoms before and after
# mutation.
_torsion_adjustments = {
    'THR': ('chi1', radians(-120)),
    'TRP': ('chi2', radians(180)),
}

refinement_starting_models = {
    'R0949': ('T0949-D1', 'T0949TS221_1-D1.pdb'),
    'R0957s2': ('T0957s2-D1', 'T0957s2TS498_1-D1.pdb'),
    'R0959': ('T0959-D1', 'T0959TS368_1-D1.pdb'),
    'R0962': ('T0962-D1', 'T0962TS246_1-D1.pdb'),
    'R0968s1': ('T0968s1-D1', 'T0968s1TS368_1-D1.pdb'),
    'R0968s2': ('T0968s2-D1', 'T0968s2TS368_3-D1.pdb'),
    'R0974s1': ('T0974s1-D1', 'T0974s1TS488_1-D1.pdb'),
    'R0976-D1': ('T0976-D1', 'T0976TS337_1-D1.pdb'),
    'R0976-D2': ('T0976-D2', 'T0976TS337_1-D2.pdb'),
    'R0977-D2': ('T0977-D2', 'T0977TS402_3-D2.pdb'),
    'R0979': ('T0979-D1', 'T0979TS470_1-D1.pdb'),
    'R0981-D3': ('T0981-D3', 'T0981TS261_1-D3.pdb'),
    'R0981-D4': ('T0981-D4','T0981TS368_1-D4.pdb'),
    'R0981-D5': ('T0981-D5', 'T0981TS116_1-D5.pdb'),
    'R0982-D2': ('T0982-D2', 'T0982TS324_1-D2.pdb'),
    'R0986s1': ('T0986s1-D1', 'T0986s1TS043_4-D1.pdb'),
    'R0986s2': ('T0986s2-D1', 'T0986s2TS043_4-D1.pdb'),
    'R0989-D1': ('T0989-D1', 'T0989TS432_2-D1.pdb'),
    'R0992': ('T0992-D1', 'T0992TS368_1-D1.pdb'),
    'R0993s2': ('T0993s2-D1', 'T0993s2TS246_1-D1.pdb'),
    'R0996-D4': ('T0996-D4', 'T0996TS324_1-D4.pdb'),
    'R0996-D5': ('T0996-D5', 'T0996TS324_1-D5.pdb'),
    'R0996-D7': ('T0996-D7', 'T0996TS324_1-D7.pdb'),
    'R0997': ('T0997-D1', 'T0997TS324_1-D1.pdb'),
    'R0999-D3': ('T0999-D3', 'T0999TS324_1-D3.pdb'),
    'R1001': ('T1001-D1', 'T1001TS368_1-D1.pdb'),
    'R1002-D2': ('T1002-D2', 'T1002TS023_1-D2.pdb'),
    'R1004-D2': ('T1004-D2', 'T1004TS324_1-D2.pdb'),
    'R1016': ('T1016-D1', 'T1016TS368_1-D1.pdb'),
}

hard_targets = ("T0954-D1", "T0957s1-D2", "T0959-D1", "T0960-D3", "T0963-D3", "T0964-D1", "T0965-D1", "T0966-D1", "T0979-D1", "T0981-D1", "T0981-D4", "T0981-D5", "T0985-D1", "T0999-D2", "T1009-D1", "T1011-D1", "T1015s2-D1", "T1021s1-D1", "T1021s2-D1", "T1022s1-D2", "T1022s2-D1")
easy_targets=("T0960-D5", "T0961-D1", "T0963-D5", "T0973-D1", "T0974s1-D1", "T0976-D1", "T0976-D2", "T0977-D1", "T0977-D2", "T0983-D1", "T0984-D1", "T0984-D2", "T0993s1-D1", "T0993s2-D1", "T0995-D1", "T0996-D1", "T0996-D2", "T0996-D3", "T0996-D4", "T0996-D5", "T0996-D6", "T0996-D7", "T0999-D3", "T0999-D4", "T0999-D5", "T1002-D1", "T1002-D2", "T1002-D3", "T1003-D1", "T1004-D1", "T1004-D2", "T1006-D1", "T1013-D1", "T1014-D1", "T1014-D2", "T1016-D1", "T1017s1-D1", "T1018-D1", "T1019s2-D1", "T1020-D1")

all_tbm_targets = hard_targets+easy_targets

tbm_gid_to_group_name = {
    322: "Zhang",
    43: "A7D",
    89: "MULTICOM",
    261: "Zhang-Server",
    460: "McGuffin",
    145: "QUARK",
    222: "Seok-refine",
    274: "MUFold",
    135: "SBROD",
    324: "RaptorX-DeepModeller",
    354: "wfAll-Cheng",
    221: "RaptorX-TBM",
    344: "Kiharalab",
    192: "Elofsson",
    197: "MESHI",
    55: "VoroMQA-select",
    86: "BAKER",
    68: "Seok",
    196: "Grudinin",
    406: "Seder3mm",
    156: "Seok-server",
    390: "Bhattacharya",
    243: "MULTICOM-CONSTRUCT",
    418: "Seder3nc",
    241: "Bhageerath-Star",
    44: "ProQ2",
    368: "BAKER-ROSETTASERVER",
    71: "Seder3full",
    441: "FALCON",
    23: "MULTICOM-NOVEL",
    366: "Venclovas",
    214: "wfRosetta-ModF7",
    208: "KIAS-Gdansk",
    58: "MULTICOM_CLUSTER",
    309: "Seder1",
    246: "IntFOLD5",
    457: "Wallner",
    426: "AP_1",
    117: "Jones-UCL",
    164: "Yang-Server",
    329: "D-Haven",
    224: "Destini",
    446: "slbio",
    116: "Zhang-CEthreader",
    281: "SHORTLE",
    471: "CPClab",
    149: "Zhou-SPOT-3D",
    266: "slbio_server",
    377: "wfRstta-Maghrabi-TQA",
    279: "ZHOU-SPOT",
    163: "Bates_BMM",
    160: "CMA-align",
    244: "Seder3hard",
    470: "Seok-assembly",
    337: "FALCON-TBM",
    335: "wfRosetta-PQ2-AngQA",
    498: "RaptorX-Contact",
    4: "YASARA",
    347: "MESHI-server",
    432: "Seok-naive_assembly",
    402: "RBO-Aleph",
    397: "PepBuilderJ",
    110: "Distill",
    47: "chuo-u",
    85: "BhageerathH-Plus",
    312: "MUFold_server",
    488: "Delta-Gelly-Server",
    112: "AWSEM",
    124: "AWSEM-Suite",
    472: "DELClab",
    358: "Spider",
    431: "Laufer",
    92: "Ricardo",
    492: "wf-BAKER-UNRES",
    288: "UNRES",
    97: "Laufer_abinitio",
    365: "3D-JIGSAW_SL1",
    407: "rawMSA",
    77: "qmo",
    351: "DL-Haven",
    157: "GAPF_LNCC",
    41: "FALCON-Contact",
    122: "Forbidden",
    7: "ACOMPMOD",
    257: "NOCONTACT",
    381: "GONGLAB-THU",
    152: "PconsC4",
    389: "UpsideUChicago",
    282: "PRAYOG",
    250: "Meilerlab",
    473: "Maurice",
    497: "GaussDCA",
    414: "BCLMeilerLab",
    401: "InnoUNRES",
    458: "FOLDNET",
    348: "HMSCasper-Refiner",
    380: "Laufer_100",
    378: "Cao-server",
    476: "Sun_Tsinghua",
}


refinement_gid_to_group_name = {
    -1:     'Starting model',
    4:      'YASARA',
    68:     'Seok',
    86:     'BAKER',
    102:    'Bhattacharya-Server',
    112:    'AWSEM',
    117:    'Jones-UCL',
    122:    'Forbidden',
    156:    'Seok-server',
    174:    'Zhang-Refinement',
    190:    'DC_refine',
    195:    'Seminoles',
    196:    'Grudinin',
    197:    'MESHI',
    208:    'KIAS-Gdansk',
    217:    'Boniecki_pred',
    270:    'Huang',
    281:    'SHORTLE',
    288:    'UNRES',
    312:    'MUFold_server',
    328:    'Kiharalab_RF2',
    329:    'D-Haven',
    344:    'Kiharalab',
    356:    'FEIGLAB',
    358:    'Spider',
    359:    '3DCNN',
    390:    'Bhattacharya',
    425:    'BAKER-AUTOREFINE',
    431:    'Laufer',
    433:    'AIR',
    457:    'Wallner',
    460:    'McGuffin',
    492:    'wf-BAKER-UNRES'
}

def get_group_names(gids, refinement=False):
    if refinement:
        name_dict = refinement_gid_to_group_name
    else:
        name_dict = tbm_gid_to_group_name
    names = [name_dict.get(gid, str(gid)) for gid in gids]
    return names



TWISTED_THRESHOLD = radians(30)

_target_dir = '/home/tic20/casp/predictioncenter.org/casp13/assessors/TARGETS'
_prediction_dir = '/home/tic20/casp/predictioncenter.org/casp13/assessors/TARBALLS/predictions_trimmed_to_domains'
_refinement_dir = '/home/tic20/casp/predictioncenter.org/casp13/assessors/TARBALLS/predictions'
_casp_analysis_dir = '/home/tic20/casp/predictioncenter.org/casp13/assessors/TARBALLS/results/tables'


def load_target(session, name):
    from chimerax.core.commands import open
    target = open.open(session, os.path.join(_target_dir, name+'.pdb'))[0]
    return target

def merge_with_casp_data(name, refinement=False, model_1_only = False):
    from collections import defaultdict
    import numpy
    db = defaultdict(list)
    model_id_to_index = dict()
    # Load the CASP standard data
    casp_filename = os.path.join(_casp_analysis_dir, name+'.txt')
    with open(casp_filename, 'rt') as f:
        data = f.read().split('\n')
    header_line = data[0].split()
    l = 0
    for line in data[2:]:
        if line:
            line = line.split()
            # Handle the first 3 columns manually
            model_id = line[1]
            if model_1_only:
                if '_1' not in model_id:
                    if model_id != 'starting_model':
                        continue
            try:
                db['#'].append(int(line[0]))
            except:
                db['#'].append(-1)
            db['Model'].append(model_id)
            model_id_to_index[model_id] = l
            group_id = line[2]
            if group_id == '-':
                db['Server'].append(False)
                group_id = -1
            elif 's' in group_id:
                db['Server'].append(True)
                group_id = int(group_id.rstrip('s'))
            else:
                db['Server'].append(False)
                group_id = int(group_id)
            db['GR#'].append(group_id)
            for h, entry in zip(header_line[3:], line[3:]):
                if entry == 'N/A' or entry =='-':
                    db[h].append(numpy.nan)
                else:
                    db[h].append(float(entry))
            l+=1
    numpy_db = dict()
    for h, d in db.items():
        numpy_db[h] = numpy.array(d)
    line_count = len(numpy_db['Model'])

    # Load the torsion-space analysis data
    if not refinement:
        base_dir = _prediction_dir
    else:
        base_dir = _refinement_dir
    torsion_filename = os.path.join(base_dir, name, name+'_torsion_comparison.csv')
    with open(torsion_filename, 'rt') as f:
        data = f.read().split('\n')
    header_line = data[0].split(',')
    for h in header_line[1:]:
        if h == 'Template Name':
            numpy_db[h] = numpy.empty(line_count, dtype=object)
        else:
            numpy_db[h] = numpy.ones(line_count)*numpy.nan
    for line in data[1:]:
        if line:
            line = line.split(',')
            model_id = line[0]
            if model_1_only:
                if '_1' not in model_id:
                    if model_id != 'starting_model':
                        continue
            try:
                index = model_id_to_index[model_id]
            except:
                print("Failed on {}".format(model_id))
                continue
            for h, entry in zip(header_line[1:], line[1:]):
                if h == 'Template Name':
                    numpy_db[h][index] = entry
                elif entry != 'N/A':
                    numpy_db[h][index] = float(entry)

    # qse = numpy_db['QSE']
    # qse[numpy.isnan(qse)] = 0
    return numpy_db

def get_value_by_name(db, name, key):
    index = numpy.argwhere(db['Model']==name)[0]
    return db[key][index]

def adjusted_zscore(data, mean=None, std=None):
    '''
    Returns an array of pseudo-z-scores processed as follows:
      - Standard Z-scores are calculated for all data
      - Data points corresponding to Z-scores less than -2 are masked out, and
        new Z-scores calculated for the remaining data. Masked-out points will
        have returned Z-scores of zero.
      - All Z-scores less than zero are set to zero.
    '''
    from scipy.stats import zscore
    z= zscore
    import numpy
    real_mask = numpy.isfinite(data)
    z_scores = numpy.zeros(len(data))
    if mean is None:
        z_scores[real_mask] = zscore(data[real_mask])
    else:
        z_scores[real_mask] = z_score_from_fixed_mean_and_std(data[real_mask], mean, std)
    z_mask = z_scores > -2
    final_mask = numpy.logical_and(z_mask, real_mask)
    if mean is None:
        z_scores[final_mask] = zscore(data[final_mask])
    else:
        z_scores[final_mask] = z_score_from_fixed_mean_and_std(data[final_mask], mean, std)
    z_scores[z_scores<0] = 0
    # z_scores[numpy.logical_not(real_mask)] = numpy.nan
    return z_scores

def z_score_from_fixed_mean_and_std(data, mean, std):
    return (data-mean)/std

def raw_score(db, names, directions, weights, ranges, mode = 'all', ranged=True):
    '''
    For each column in names, the scores are divided by their range, clamped
    to (0..1). If the direction is -1, the result is subtracted from 1.
    The resulting values are averaged using the given weights.
    '''
    import numpy
    scores = []
    coverage = db['Model Coverage']
    gids = db['GR#']
    _, unique_indices = numpy.unique(gids, return_index=True)
    unique_gids = gids[numpy.sort(unique_indices)]
    for name, d, w, r in zip(names, directions, weights, ranges):
        data = db[name].copy()
        if ranged:
            data /= r
            data[data<0] = 0
            data[data>1] = 1
            if d == -1:
                data = 1-data
        elif d == -1:
            data = -data
        data[numpy.logical_not(numpy.isfinite(data))]=0
        scores.append(data*w)
    result = sum(scores)/sum(weights)
    if mode == 'all':
        return gids, result*coverage
    elif mode == 'model 1 only':
        model_ids = db['Model']
        model_mask = numpy.array([True if '_1' in mid else False for mid in model_ids])
        model_mask[db['GR#']==-1] = True
        return db['GR#'][model_mask], result[model_mask] * coverage[model_mask]
    elif mode == 'best only':
        final_result = numpy.empty(len(unique_gids))
        for i, gid in enumerate(unique_gids):
            gmask = (gids==gid)
            final_result[i] = numpy.max(result[gmask]*coverage[gmask])
        return unique_gids, final_result

def weighted_z_score(db, names, directions, weights, adjusted=True,
        mode='all'):
    '''
    For each column in names, calculates a Z-score for each point based on the
    mean and standard deviation for that column. Combines the Z-scores according
    to the given weights, and optionally adjusts them according to CASP standard
    procedure using adjusted_zscore(). 'directions' should be a list containing
    strictly +/-1, where +1 indicates a higher score for the given criterion is
    better, and -1 indicates a lower score is better
    '''
    import numpy
    z_scores = []
    coverage = db['Model Coverage']
    gids = db['GR#']
    gid_mask = (gids != -1)
    _, unique_indices = numpy.unique(gids, return_index=True)
    unique_gids = gids[numpy.sort(unique_indices)]

    for name, d, w in zip(names, directions, weights):
        data = db[name]*d
        fmask = numpy.logical_and(gid_mask, numpy.isfinite(data))
        mean = numpy.mean(data[fmask])
        std = numpy.std(data[fmask])
        if adjusted:
            # z_score = adjusted_zscore(data, mean=mean, std=std)
            z_score = adjusted_zscore(data)
        else:
            from scipy.stats import zscore
            real_mask = numpy.isfinite(data)
            fmask = numpy.logical_and(real_mask, gid_mask)
            z_score = numpy.zeros(len(data))
            z_score[fmask] = zscore(data[fmask])
        z_scores.append(z_score*w)
    result = sum(z_scores) / sum(weights)
    if mode in ('all', 'model_1_only'):
        return gids, result*coverage

    final_result = numpy.empty(len(unique_gids))
    # if mode == 'model 1 only':
    #     model_ids = db['Model']
    #     model_mask = numpy.array([True if '_1' in mid else False for mid in model_ids ])
    #     model_mask[db['GR#']==-1] = True
    #     final_result = numpy.empty(sum(model_mask))
    #     try:
    #         final_result[:] = result[model_mask] * coverage[model_mask]
    #     except:
    #         print(db['Model'][0], numpy.sort(unique_gids), numpy.sort(db['GR#'][model_mask]))
    #         raise
    #     return db['GR#'][model_mask], final_result
    if mode == 'best only':
        for i, gid in enumerate(unique_gids):
            gmask = (gids==gid)
            final_result[i] = numpy.max(result[gmask] * coverage[gmask])
    final_result[numpy.isnan(final_result)] = 0
    return unique_gids, final_result

def zero():
    return 0

def choose_mode(model_1_only, best_only):
    if model_1_only and best_only:
        raise TypeError('Cannot select both model_1_only and best_only!')
    if model_1_only:
        return 'model 1 only'
    elif best_only:
        return 'best only'
    else:
        return 'all'

def rank_scores(group_ids, scores, plot=True):
    '''
    scores should be in the same order as the data in the db. For each group,
    the scores from submitted models will be averaged, and optionally plotted as
    a bar chart against group id in decreasing order of score.
    '''
    from collections import defaultdict
    import numpy
    sums = defaultdict(zero)
    counts = defaultdict(zero)
    for gid, score in zip(group_ids, scores):
        sums[gid]+=score
        counts[gid]+=1
    gids = numpy.array(list(sums.keys()), int)
    g_scores = numpy.empty(len(gids))
    for i, id in enumerate(gids):
        g_scores[i] = sums[id]/counts[id]
    order = numpy.argsort(g_scores)[::-1]
    gids[:] = gids[order]
    g_scores[:] = g_scores[order]

    if plot:
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x = numpy.arange(len(gids))
        ax.bar(x, g_scores)
        ax.set_xticks(x)
        ax.set_xticklabels(gids, rotation=45, fontsize='small')
        plt.show()
    return gids, g_scores

def _scoring_base(db, names, directions, weights, ranges, z_scores=True, mode='all', ranged=True):
    if z_scores:
        gids, scores = weighted_z_score(db, names, directions, weights, mode=mode )
    else:
        gids, scores = raw_score(db, names, directions, weights, ranges, mode=mode, ranged=ranged)
    return gids, scores

def single_criterion_fn(label, direction=1, max_val=1):
    def f(db, z_scores=True, mode='all', plot=True):
        gids, scores = _scoring_base(db, [label], [direction], [1], [max_val],
            z_scores=z_scores, mode=mode, ranged=False)
        return gids, scores
        # return rank_scores(gids, scores, plot)
    if direction < 0:
        f.__name__='minus_{}'.format(label)
    else:
        f.__name__=label
    return f

def default_tbm_ranking(db, z_scores=True, mode='all', plot=True):
    names=['GDT_HA', 'LDDT', 'CAD_AA', 'SphGr', 'QSE']
    directions=[1,1,1,1,1]
    weights=[1/3,1/9,1/9,1/9,1/3]
    ranges=[100, 1, 1, 100, 100]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           z_scores=z_scores, mode=mode)
    return gids, scores
    #return rank_scores(gids, scores, plot)

def default_tbm_minus_qse_ranking(db, z_scores=True, mode='all', plot=True):
    names=['GDT_HA', 'LDDT', 'CAD_AA', 'SphGr']
    directions=[1,1,1,1]
    weights=[1/3,1/9,1/9,1/9]
    ranges=[100, 1, 1, 100]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           z_scores=z_scores, mode=mode)
    return rank_scores(gids, scores, plot)


def default_refinement_ranking(db, z_scores=True, mode='all', plot=True):
    names = ['RMS_CA', 'GDT_HA', 'SphGr', 'QCS', 'MolPrb_Score']
    directions = [-1, 1, 1, 1, -1]
    weights = [0.46, 0.17, 0.2, 0.15, 0.02]
    ranges=[10, 100, 100, 100, 100, 5]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           z_scores=z_scores, mode=mode)
    return rank_scores(gids, scores, plot)

def geom_quality_ranking(db, z_scores=True, mode='all', plot=True):
    names=['GDT_HA','LDDT', 'CAD_AA', 'SphGr', 'MolPrb_clash', 'Model Average Backbone Error', 'Model Average Sidechain Error', 'QSE']
    directions=[1,1,1,1,-1,-1,-1,1]
    weights=[1/4,1/16,1/16,1/16,1/8,1/8,1/16,1/4]
    ranges=[100, 1, 1, 100, 20, 45, 45, 100]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           z_scores=z_scores, mode=mode)
    return rank_scores(gids, scores, plot)

def geom_quality_ranking_new(db, z_scores=True, mode='all', plot=True):
    names=['GDT_HA','LDDT', 'CAD_AA', 'SphGr', 'MolPrb_clash', 'Model Backbone Score', 'Model Sidechain Score', 'QSE']
    directions=[1,1,1,1,-1,-1,-1,1]
    weights=[1/4,1/16,1/16,1/16,1/8,1/8,1/16,1/4]
    ranges=[100, 1, 1, 100, 20, 1, 1, 100]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           z_scores=z_scores, mode=mode)
    return rank_scores(gids, scores, plot)

def geom_quality__minus_qse_ranking(db, z_scores=True, mode='all', plot=True):
    names=['GDT_HA','LDDT', 'CAD_AA', 'SphGr', 'MolPrb_clash', 'Model Average Backbone Error', 'Model Average Sidechain Error']
    directions=[1,1,1,1,-1,-1,-1]
    weights=[1/3,1/12,1/12,1/12,1/6,1/6,1/12]
    ranges=[100, 1, 1, 100, 20, 45, 45]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           z_scores=z_scores, mode=mode)
    return rank_scores(gids, scores, plot)

def geom_quality__minus_qse_ranking_new(db, z_scores=True, mode='all', plot=True):
    names=['GDT_HA','LDDT', 'CAD_AA', 'SphGr', 'MolPrb_clash', 'Model Backbone Score', 'Model Sidechain Score']
    directions=[1,1,1,1,-1,-1,-1]
    weights=[1/3,1/12,1/12,1/12,1/6,1/6,1/12]
    ranges=[100, 1, 1, 100, 20, 1, 1]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           z_scores=z_scores, mode=mode)
    return rank_scores(gids, scores, plot)

def pure_torsion_ranking(db, z_scores=True, mode='all', plot=True):
    names=['Model Average Backbone Error', 'Model Average Sidechain Error']
    directions=[-1,-1]
    weights=[2/3,1/3]
    ranges=[45, 45]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           mode=mode, z_scores=z_scores)
    return rank_scores(gids, scores, plot)

def pure_torsion_ranking_new(db, z_scores=True, mode='all', plot=True):
    names=['Model Backbone Score', 'Model Sidechain Score']
    directions=[-1,-1]
    weights=[2/3,1/3]
    ranges=[1, 1]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                           mode=mode, z_scores=z_scores)
    return rank_scores(gids, scores, plot)


def torsion_plus_clash_ranking(db, z_scores=True, mode='all', plot=True):
    names=['Model Average Backbone Error', 'Model Average Sidechain Error', 'MolPrb_clash']
    directions=[-1,-1,-1]
    weights=[1,1,1]
    ranges = [45, 45, 20]
    gids, scores = _scoring_base(db, names, directions, weights, ranges,
                                 mode=mode, z_scores=z_scores)
    return rank_scores(gids, scores, plot)

def rmsd_ca(db, z_scores=True, mode='all', plot=True):
    if z_scores:
        gids, scores = _scoring_base(db, ['RMS_CA'], [-1], [1], [10],
                                     mode=mode, z_scores=True)
    else:
        scores = db['RMS_CA'].copy()
        gids = db['GR#'].copy()
        import numpy
        _, unique_indices = numpy.unique(gids, return_index=True)
        unique_gids = gids[numpy.sort(unique_indices)]
        result = numpy.empty(len(unique_gids))
        if mode == 'model 1 only':
            mmask = numpy.array([True if '_1' in g else False for g in db['Model'] ])
            mmask[gids==-1] = True
            result[:] = scores[mmask]
            gids, scores = unique_gids, result
        elif mode == 'best only':
            for i, gid in enumerate(unique_gids):
                result[i] = max(scores[gids==gid])
            gids, scores = unique_gids, result
        elif mode != 'all':
            raise TypeError('Unknown scoring mode!')
    return rank_scores(gids, scores, plot)

def refinement_improvements_over_starting_model(targets, scoring_fn,
        model_1_only=False, best_only=False, plot=True, z_scores=True,
        group_names=True, num_top_groups=None):
    from collections import defaultdict
    import numpy
    gid_to_num_improved = defaultdict(zero)
    gid_to_scores = defaultdict(list)
    for t in targets:
        gids, scores = refinement_improvement(t, scoring_fn,
            model_1_only=model_1_only, best_only=best_only, plot=False,
            z_scores=z_scores)
        for gid, score in zip(gids, scores):
            if score > 0:
                gid_to_num_improved[gid] += 1
            if numpy.isfinite(score):
                gid_to_scores[gid].append(score)
    final_gids = numpy.array(list(gid_to_scores.keys()))
    final_scores = numpy.zeros(len(final_gids))
    for i, gid in enumerate(final_gids):
        scores = gid_to_scores[gid]
        if len(scores) == 0:
            avg = 0
        else:
            avg = sum(scores)
        # if gid in (68, 156):
        #     print('Output {}: {}'.format(gid, avg)) #','.join('{:.2f}'.format(score) for score in scores)))
        final_scores[i] = avg

    final_counts = []
    for gid in final_gids:
        final_counts.append(gid_to_num_improved[gid])
    final_counts = numpy.array(final_counts)

    final_score_order = numpy.argsort(final_scores)[::-1]
    final_count_order = numpy.argsort(final_counts)[::-1]

    total_z_scores = (final_gids[final_score_order], final_scores[final_score_order])
    improved_counts = (final_gids[final_count_order], final_counts[final_count_order])

    if not num_top_groups:
        num_top_groups = len(final_gids)

    if plot:
        from matplotlib import pyplot as plt
        fig, axs = plt.subplots(2)
        # fig.suptitle('Refinement improvements over starting model by {}'.format(scoring_fn.__name__))
        x = numpy.arange(num_top_groups)
        ax = axs[0]
        ax.bar(x, total_z_scores[1][:num_top_groups])
        ax.set_xticks(x)
        if group_names:
            labels = get_group_names(total_z_scores[0][:num_top_groups], refinement=True)
        else:
            labels = total_z_scores[0][:num_top_groups]
        ax.set_xticklabels(labels, rotation=90)#, fontsize='small')
        ax.set_ylabel('Σ(ΔZ)')

        ax = axs[1]
        ax.bar(x, improved_counts[1][:num_top_groups]/len(targets)*100)
        ax.set_xticks(x)
        if group_names:
            labels = get_group_names(improved_counts[0][:num_top_groups], refinement=True)
        else:
            labels = improved_counts[0][:num_top_groups]
        ax.set_xticklabels(labels, rotation=90)#, fontsize='small')
        ax.set_ylabel('Improved models (%)')
        plt.tight_layout()
        plt.show()
        return plt
    return (total_z_scores, improved_counts)


def refinement_improvement(target, scoring_fn, model_1_only=False, best_only=False, plot=True,
        z_scores=False):
    mode = choose_mode(model_1_only, best_only)
    import numpy
    db = merge_with_casp_data(target, refinement=True, model_1_only=model_1_only)
    gids, scores = scoring_fn(db, z_scores=z_scores, mode=mode, plot=False)
    # print('Reference score: {}'.format(scores[gids==-1][0]))
    # print('Input scores before correction: 68: {:.2f}; 156: {:.2f}'.format(scores[gids==68][0], scores[gids==156][0]))
    scores -= scores[gids==-1][0]
    gid_mask = gids!=-1
    gids=gids[gid_mask]
    scores = scores[gid_mask]
    # print('Input scores: 68: {:.2f}; 156: {:.2f}'.format(scores[gids==68][0], scores[gids==156][0]))
    if plot:
        from matplotlib import pyplot as plt
        fig = plt.figure()
        fig.suptitle('{} refinement improvements over starting model ({})'.format(target, mode))
        ax = fig.add_subplot(111)
        x = numpy.arange(len(scores))
        ax.bar(x, scores)
        ax.set_xticks(x)
        ax.set_xticklabels(gids, rotation=45, fontsize='small')
        ax.set_ylabel('Raw {} score'.format(scoring_fn.__name__))
        plt.show()
    return gids, scores

def refinement_improvement_histogram(targets, scoring_fn, model_1_only=False, best_only=False,
        n_bins=20, z_scores=True, range = (-10,10)):
    mode = choose_mode(model_1_only, best_only)
    import numpy
    collated_scores = []
    for t in targets:
        _, scores = refinement_improvement(t, scoring_fn, model_1_only=model_1_only,
                best_only=best_only, z_scores=z_scores, plot=False)
        collated_scores.extend(list(scores))
    frac_above_zero = sum([True for s in collated_scores if s > 0]) / len(collated_scores)
    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    final_scores = []
    for s in collated_scores:
        if s>range[0] and s<range[1]:
            final_scores.append(s)
    ax.hist(final_scores, bins=n_bins)
    ax.set_xlabel('Δ({})'.format(scoring_fn.__name__))
    ax.set_ylabel('Frequency')
    ax.text(0.025, 0.9, 'Improved over template: {:.2f}%'.format(frac_above_zero*100),
        horizontalalignment='left', transform=ax.transAxes)
    plt.show()


def compare_refinement_improvements(target, scoring_fn_1, scoring_fn_2,
        model_1_only=False, best_only=False, plot=True):
    import numpy
    x_gids, x_scores = refinement_improvement(target, scoring_fn_1,
        model_1_only=model_1_only, best_only=best_only, plot=False)
    y_gids, y_scores = refinement_improvement(target, scoring_fn_2,
        model_1_only=model_1_only, best_only=best_only, plot=False)

    shifts = ranking_shift(y_gids, x_gids)
    order_x = numpy.argsort(x_gids)
    gids = x_gids[order_x]
    x_scores=x_scores[order_x]

    order_y = numpy.argsort(y_gids)
    y_scores = y_scores[order_y]
    shifts = shifts[order_y]

    nan_mask = numpy.logical_and(numpy.isfinite(x_scores), numpy.isfinite(y_scores))
    x_scores = x_scores[nan_mask]
    x_gids = x_gids[nan_mask]
    y_scores = y_scores[nan_mask]
    y_gids = y_gids[nan_mask]
    shifts = shifts[nan_mask]

    if plot:
        ranking_comparison_plot(x_scores, y_scores, shifts, gids, scoring_fn_1.__name__,
            scoring_fn_2.__name__, refinement=True)
    return (x_scores, y_scores, shifts, gids)

def _ninetieth_percentile(data):
    import numpy
    return numpy.percentile(data, 90)

def refinement_improvement_by_model_source(targets, scoring_fn, model_1_only=False,
        best_only=False, plot=True, mode='best'):
    import numpy
    if mode == 'best':
        stat_fn = numpy.max
    elif mode == 'mean':
        stat_fn = numpy.mean
    elif mode == 'median':
        stat_fn = numpy.median
    elif mode == '90th percentile':
        stat_fn = _ninetieth_percentile
    else:
        raise TypeError('mode argument must be one of "best", "mean" or "median"')
    from collections import defaultdict
    score_dict = defaultdict(zero)
    counts = defaultdict(zero)
    for t in targets:
        template_name = refinement_starting_models[t][1]
        template_gid = template_name.split('_')[0][-3:]
        *_, scores = refinement_improvement(t, scoring_fn,
            model_1_only=model_1_only, best_only=best_only, plot=False)
        nan_mask = numpy.isfinite(scores)
        score_dict[template_gid] += stat_fn(scores[nan_mask])
        counts[template_gid] += 1
    gids = list(score_dict.keys())
    final_scores = []
    for gid in gids:
        final_scores.append(score_dict[gid]/counts[gid])

    if plot:
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x = numpy.arange(len(gids))
        ax.bar(x, final_scores)
        ax.set_xticks(x)
        ax.set_xticklabels(gids, rotation=45, fontsize='small')
        ax.set_xlabel('Source of starting model')
        ax.set_ylabel('Average of {} raw {} values'.format(mode, scoring_fn.__name__))
        plt.show()
    return gids, final_scores

def refinement_improvement_by_starting_quality(targets, x_scoring_fn, y_scoring_fn, model_1_only=False,
    best_only=False, plot=True, mode='best', reverse=False):
    scoring_mode = choose_mode(model_1_only, best_only)
    import numpy
    if mode == 'best':
        stat_fn = numpy.max
    elif mode == 'mean':
        stat_fn = numpy.mean
    elif mode == 'median':
        stat_fn = numpy.median
    elif mode == '90th percentile':
        stat_fn = _ninetieth_percentile
    else:
        raise TypeError('mode argument must be one of "best", "mean" or "median"')
    starting_scores = []
    final_scores=[]
    for t in targets:
        db = merge_with_casp_data(t, refinement=True)
        gids, raw_scores = x_scoring_fn(db, z_scores=False, mode=scoring_mode, plot=False)
        starting_index = numpy.argwhere(gids==-1)[0][0]
        starting_score = raw_scores[starting_index]
        starting_scores.append(starting_score)
        *_, improvements = refinement_improvement(t, y_scoring_fn,
            model_1_only=model_1_only, best_only=best_only,
            plot=False)
        if not reverse:
            final_scores.append(stat_fn(improvements))
        else:
            final_scores.append(-stat_fn(-improvements))

    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(starting_scores, final_scores)


    if plot:
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(starting_scores, final_scores)
        ax.set_xlabel('Starting raw {} score'.format(x_scoring_fn.__name__))
        ax.set_ylabel('{} improvement raw {} score'.format(mode, y_scoring_fn.__name__))
        lr_xvals = numpy.linspace(min(starting_scores), max(starting_scores), num=2)
        lr_yvals = slope*lr_xvals + intercept
        fig.suptitle('{} refinement improvement vs. starting model score (R**2 = {:0.2f})'.format(mode, r_value**2))
        ax.plot(lr_xvals, lr_yvals)

        plt.show()
    return starting_scores, final_scores

def refinement_improvement_by_secondary_structure_content(session, targets, scoring_fn,
        model_1_only=False, best_only=False, plot=True, mode='best'):
    import numpy
    if mode == 'best':
        stat_fn = numpy.max
    elif mode == 'mean':
        stat_fn = numpy.mean
    elif mode == 'median':
        stat_fn = numpy.median
    elif mode == '90th percentile':
        stat_fn = _ninetieth_percentile
    else:
        raise TypeError('mode argument must be one of "best", "mean" or "median"')
    xvals = []
    yvals = []
    for t in targets:
        tm = load_target(session, t)
        from chimerax.atomic import Residue
        res = tm.residues
        ss_frac = sum(res.ss_types != Residue.SS_COIL)/len(res)
        xvals.append(ss_frac)
        session.models.close([tm])
        *_, improvements = refinement_improvement(t, scoring_fn,
            model_1_only=model_1_only, best_only=best_only,
            plot=False)
        yvals.append(stat_fn(improvements))

    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(xvals, yvals)

    if plot:
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(xvals, yvals)
        lr_xvals = numpy.linspace(min(xvals), max(xvals), num=2)
        lr_yvals = slope*lr_xvals + intercept
        fig.suptitle('{} refinement improvement vs. 2ry structure content (R**2 = {:0.2f})'.format(mode, r_value**2))
        ax.plot(lr_xvals, lr_yvals)
        ax.set_xlabel('Helix or beta sheet content')
        ax.set_ylabel('{} improvement in raw {} score'.format(mode, scoring_fn.__name__))
        plt.show()
    return xvals, yvals, r_value**2

def refinement_before_vs_after_scores(targets, scoring_fn, num_top_groups=5,
        target_groups=None, reverse=False, model_1_only=False,
        best_only=False, plot=True):
    mode = choose_mode(model_1_only, best_only)
    import numpy
    if target_groups is None:
        target_dict = {}
        from collections import defaultdict
        all_gids=set()
        for t in targets:
            gid_to_scores = target_dict[t] = defaultdict(list)
            db = merge_with_casp_data(t, refinement=True)
            gids, scores = scoring_fn(db, z_scores=True, mode=mode, plot=False)
            for gid, score in zip(gids, scores):
                if gid != -1:
                    gid_to_scores[gid].append(score)
                    all_gids.add(gid)
        z_score_dict = {}
        for gid in all_gids:
            z_scores = []
            for g2s in target_dict.values():
                g_scores = g2s[gid]
                if len(g_scores) < 5:
                    for i in range(5-len(g_scores)):
                        g_scores.append(0)
                z_scores.extend(g_scores)
            z_score_dict[gid] = numpy.array(z_scores).mean()
        gids = []
        z_scores = []
        for gid, z_score in z_score_dict.items():
            gids.append(gid)
            z_scores.append(z_score)
        sorted_gids = [gid for gid, _ in sorted(zip(gids, z_scores), key=lambda pair: pair[1], reverse=not reverse)]
        target_groups = sorted_gids[:num_top_groups]
    final_score_dict = {}
    for gid in target_groups:
        final_score_dict[gid] = { 'x': [], 'y': [], 'target': []}
    for t in targets:
        db = merge_with_casp_data(t, refinement=True)
        gids, scores = scoring_fn(db, z_scores=False, mode=mode, plot=False)
        starting_index = numpy.argwhere(gids==-1).ravel()[0]
        starting_score = scores[starting_index]
        for gid in target_groups:
            final_indices = numpy.argwhere(gids==gid).ravel()
            for s in scores[final_indices]:
                final_score_dict[gid]['x'].append(starting_score)
                final_score_dict[gid]['y'].append(s)
                final_score_dict[gid]['target'].append(t)


    if plot:
        from matplotlib import pyplot as plt
        from matplotlib import cm
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.suptitle('Initial vs. final {} for top {} groups'.format(scoring_fn.__name__, num_top_groups))

        colors = cm.rainbow(numpy.linspace(0, 1, num_top_groups))
        max_score = max([max([max(vals) for vals in (d['x'], d['y'])]) for d in final_score_dict.values()])
        min_score = min([min([min(vals) for vals in (d['x'], d['y'])]) for d in final_score_dict.values()])
        scatter_plots = {}
        for (gid, data), c in zip(final_score_dict.items(), colors):
            scatter_plots[gid] = ax.scatter(data['x'], data['y'], color=c, label='Group {}'.format(gid))
        ax.legend()
        ax.plot((min_score, max_score), (min_score, max_score))
        ax.set_xlabel('Starting {} score'.format(scoring_fn.__name__))
        ax.set_ylabel('Final {} score'.format(scoring_fn.__name__))

        annot = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords='offset points',
                    bbox=dict(boxstyle='round', fc='w'),
                    arrowprops=dict(arrowstyle='->'))
        annot.set_visible(False)

        def update_annot(ind_dict, pos):
            text = ''
            for gid, indices in ind_dict.items():
                targets = (final_score_dict[gid]['target'][i] for i in indices)
                text += '{}: {}\n'.format(gid, ', '.join(targets))
            annot.set_text(text)
            annot.xy = pos

        def hover(event):
            if event.inaxes == ax:
                ind_dict = {}
                for gid, sc in scatter_plots.items():
                    cont, ind = sc.contains(event)
                    if cont:
                        indices = ind['ind']
                        ind_dict[gid] = indices
                        pos = sc.get_offsets()[indices][0]
                if ind_dict:
                    update_annot(ind_dict, pos)
                    annot.set_visible(True)
                else:
                    annot.set_visible(False)
                fig.canvas.draw()

        fig.canvas.mpl_connect('motion_notify_event', hover)

        plt.show()

    return final_score_dict






def ranking_dict(targets, scoring_fn, model_1_only=False, best_only=False, refinement=False):
    mode = choose_mode(model_1_only, best_only)
    target_rankings = dict()
    for t in targets:
        db = merge_with_casp_data(t, refinement=refinement, model_1_only=model_1_only)
        try:
            target_rankings[t] = scoring_fn(db, mode=mode, plot=False)
        except KeyError:
            raise KeyError('Error on target {}'.format(t))
    return target_rankings

def meta_ranking(targets, scoring_fn, plot=True, num_groups=None,
    model_1_only=False, best_only=False, refinement=False, use_group_names=False,
    rotation=45):
    '''
    Rank groups by their average score for all targets in the given list, using
    the chosen scoring function. Set num_groups to a positive value to plot only
    that number of leaders.
    '''
    import os
    import numpy
    target_rankings = ranking_dict(targets, scoring_fn, model_1_only=model_1_only,
        best_only=best_only, refinement=refinement)

    from collections import defaultdict
    sums=defaultdict(zero)
    count=0
    for t, data in target_rankings.items():
        for gid, g_score in zip(*data):
            if not numpy.isnan(g_score):
                sums[gid] += g_score
        count += 1
    gids = numpy.array(list(sums.keys()), int)
    g_scores = numpy.empty(len(gids))
    for i, id in enumerate(gids):
        g_scores[i]=sums[id] # /count
    order = numpy.argsort(g_scores)[::-1]
    gids[:] = gids[order]
    g_scores[:]=g_scores[order]

    p_gids = gids
    p_gscores = g_scores
    if num_groups is not None:
        p_gids = p_gids[:num_groups]
        p_gscores = p_gscores[:num_groups]
    if plot:
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x = numpy.arange(len(p_gids))
        if refinement:
            colors = numpy.empty([len(x), 3])
            colors[:] = numpy.array([31, 119, 180])/255
            colors[p_gids == -1] = [1.0,0,0]
            ax.bar(x, p_gscores, color=colors)
        else:
            ax.bar(x, p_gscores)
        if not use_group_names:
            labels = p_gids
        else:
            labels = get_group_names(p_gids, refinement=refinement)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=rotation, fontsize='small')
        ax.set_ylabel(scoring_fn.__name__)
        plt.tight_layout()
        plt.show()
        return plt
    return gids, g_scores

def per_model_rankings(gid, targets, scoring_fns, model_1_only=False,
        best_only=False, refinement=False, plot=True):
    '''
    Gives the ranking of the chosen group for each target in targets, using the
    given ranking function.
    '''
    import numpy
    rankings_by_function = {}
    targets = numpy.array(targets)
    for scoring_fn in scoring_fns:
        target_rankings = ranking_dict(targets, scoring_fn, model_1_only=model_1_only,
            best_only=best_only,
            refinement=refinement)
        z_scores = numpy.ones(len(targets))*numpy.nan
        for t, data in target_rankings.items():
            gids, g_scores = data
            indices = numpy.argwhere(gids==gid)
            if len(indices):
                z_scores[targets==t] = numpy.mean(g_scores[indices])
        rankings_by_function[scoring_fn.__name__] = z_scores
    if plot:
        from matplotlib import pyplot as plt
        fig, axs = plt.subplots(len(rankings_by_function), sharex=True)
        # If only one scoring function is given, axs is a AxesSubplot instance
        # rather than a list, necessitating the below hack
        if not hasattr(axs, '__len__'):
            axs=[axs]
        fig.suptitle('Per-model rankings for group {}'.format(gid))
        x = numpy.arange(len(targets))
        for i, (fn, scores) in enumerate(rankings_by_function.items()):
            ax = axs[i]
            nan_mask = numpy.isnan(scores)
            colors = numpy.empty([len(x), 3])
            colors[:] = numpy.array([31, 119, 180])/255
            colors[nan_mask] = numpy.array([220, 220, 220])/255
            scores[nan_mask] = 1.0
            ax.bar(x, scores, color=colors)
            ax.set_xticks(x)
            ax.set_ylim([0,1.5])
            ax.set_title(fn)
        ax.set_xticklabels(targets, rotation=45, fontsize='small')
        plt.show()
    return rankings_by_function






def ranking_shift(sorted_gids_1, sorted_gids_2):
    '''
    Returns an array giving the number of ranks each id in sorted_gids_1 has
    improved/regressed by in sorted_gids_2
    '''
    import numpy
    def invperm(p):
        q = numpy.empty_like(p)
        q[p] = numpy.arange(len(p))
        return q

    o1 = numpy.argsort(sorted_gids_1)
    o2 = numpy.argsort(sorted_gids_2)
    shifts = numpy.arange(len(o1)) - o2[invperm(o1)]
    return shifts


def compare_rankings(db, rank_fn_1, rank_fn_2, model_1_only=False, best_only=False,
        plot=True, refinement=False, use_group_names=True):
    mode = choose_mode(model_1_only, best_only)
    d_gids, d_scores = rank_fn_1(db, mode=mode, plot=False)
    g_gids, g_scores = rank_fn_2(db, mode=mode, plot=False)
    import numpy
    orderd = numpy.argsort(d_gids)
    sorted_gids = d_gids[orderd]
    sorted_shifts = shifts[orderd]
    sorted_d_scores = d_scores[orderd]
    orderg = numpy.argsort(g_gids)
    sorted_g_scores = g_scores[orderg]
    if plot:
        ranking_comparison_plot(sorted_d_scores, sorted_g_scores, -sorted_shifts,
            sorted_gids, rank_fn_1.__name__, rank_fn_2.__name__, refinement=refinement, use_group_names=use_group_names)


def compare_meta_rankings(targets, rank_fn_1, rank_fn_2, num_groups=None,
    model_1_only = False, best_only=False, refinement=False, use_group_names=True):
    x_gids, x_scores = meta_ranking(targets, rank_fn_1, plot=False,
        model_1_only=model_1_only, best_only=best_only, refinement=refinement)
    y_gids, y_scores = meta_ranking(targets, rank_fn_2, plot=False,
        model_1_only=model_1_only, best_only=best_only, refinement=refinement)

    shifts = ranking_shift(y_gids, x_gids)
    import numpy
    order_x = numpy.argsort(x_gids)
    gids = x_gids[order_x]
    x_scores = x_scores[order_x]

    order_y = numpy.argsort(y_gids)
    y_scores = y_scores[order_y]
    shifts = shifts[order_y]

    nan_mask = numpy.logical_not(numpy.logical_or(numpy.isnan(x_scores), numpy.isnan(y_scores)))
    x_scores = x_scores[nan_mask]
    x_gids = x_gids[nan_mask]
    y_scores = y_scores[nan_mask]
    y_gids = y_gids[nan_mask]
    shifts = shifts[nan_mask]

    if num_groups is not None:
        gids = gids[:num_groups]
        x_scores = x_scores[:num_groups]
        y_scores = y_scores[:num_groups]
        shifts = shifts[:num_groups]

    ranking_comparison_plot(x_scores, y_scores, shifts, gids, rank_fn_1.__name__, rank_fn_2.__name__, refinement=refinement, use_group_names=use_group_names)
    return (x_scores, y_scores, shifts, gids)

def ranking_comparison_plot(xscores, yscores, shifts, gids, x_fn_name, y_fn_name, refinement=False, use_group_names=True):
    import numpy
    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(x_fn_name)
    ax.set_ylabel(y_fn_name)

    sp = ax.scatter(xscores, yscores, c=shifts, cmap='cool')

    if refinement:
        s_mask = gids==-1
        ax.scatter(xscores[s_mask], yscores[s_mask], c='red', s=49)

    if use_group_names:
        if refinement:
            name_dict = refinement_gid_to_group_name
        else:
            name_dict = tbm_gid_to_group_name

    annot = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords='offset points',
                bbox=dict(boxstyle='round', fc='w'),
                arrowprops=dict(arrowstyle='->'))
    annot.set_visible(False)

    def update_annot(indices, xdata, ydata):
        xs = xscores[indices]
        ys = yscores[indices]
        nan_mask = numpy.logical_not(numpy.logical_or(numpy.isnan(xs), numpy.isnan(ys)))
        xs = xs[nan_mask]
        ys = ys[nan_mask]
        if not len(xs):
            return False
        x = numpy.mean(xs)
        y = numpy.mean(ys)
        lgids = gids[indices]
        if use_group_names:
            names = get_group_names(lgids, refinement=refinement)
        else:
            names = lgids
        text = ', '.join([str(name) for name in names])
        annot.set_text(text)
        annot.xy = x, y
        return True

    def hover(event):
        if event.inaxes == ax:
            cont, ind = sp.contains(event)
            if cont:
                candidates = ind['ind']
                result = update_annot(candidates, xscores, yscores)
                annot.set_visible(result)
            else:
                annot.set_visible(False)
            fig.canvas.draw()

    fig.canvas.mpl_connect('motion_notify_event', hover)

    plt.show()


# def plot_geom_ranking_vs_default(db):
#     d_gids, d_scores = default_ranking(db, plot=False)
#     g_gids, g_scores = geom_quality_ranking(db, plot=False)
#     shifts = ranking_shift(d_gids, g_gids)
#     import numpy
#     orderd = numpy.argsort(d_gids)
#     sorted_gids = d_gids[orderd]
#     sorted_shifts = shifts[orderd]
#     sorted_d_scores = d_scores[orderd]
#     orderg = numpy.argsort(g_gids)
#     sorted_g_scores = g_scores[orderg]
#
#
#     from matplotlib import pyplot as plt
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.set_xlabel('Default ranking score')
#     ax.set_ylabel('Local conformation score')
#     sp = ax.scatter(sorted_d_scores, sorted_g_scores, c=-sorted_shifts, cmap='cool')
#
#     annot = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords='offset points',
#                 bbox=dict(boxstyle='round', fc='w'),
#                 arrowprops=dict(arrowstyle='->'))
#     annot.set_visible(False)
#
#     def update_annot(indices, xdata, ydata):
#         x = numpy.mean(sorted_d_scores[indices])
#         y = numpy.mean(sorted_g_scores[indices])
#         names = sorted_gids[indices]
#         text = ', '.join([str(name) for name in names])
#         annot.set_text(text)
#         annot.xy = x, y
#
#     def hover(event):
#         if event.inaxes == ax:
#             cont, ind = sp.contains(event)
#             if cont:
#                 candidates = ind['ind']
#                 xdata = sorted_d_scores
#                 ydata = sorted_g_scores
#                 update_annot(candidates, xdata, ydata)
#                 annot.set_visible(True)
#             else:
#                 annot.set_visible(False)
#             fig.canvas.draw()
#
#     fig.canvas.mpl_connect('motion_notify_event', hover)
#
#     plt.show()



def _xyz_from_string(xyz):
    '''
    Ugly hack to get xyz coordinate from the string generated by Axes3D's
    format_coord() function.
    '''
    xyz_string = xyz.split(',')
    xyz_float = [float(s[2:]) for s in xyz_string]
    return xyz_float


def plot_3d(db, labels, server_only = False, best_only = False, limits=[None, None, None]):
    from mpl_toolkits.mplot3d import Axes3D, proj3d
    from matplotlib import pyplot as plt
    import numpy
    s_mask = db['Server'].copy()
    if not server_only:
        s_mask[:] = True
    b_mask = numpy.ones(len(s_mask), numpy.bool)
    if best_only:
        names = db['Model']
        for i, name in enumerate(names):
            if '_1' not in name:
                b_mask[i] = False
    final_mask = numpy.logical_and(s_mask, b_mask)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sp = ax.scatter(db[labels[0]][final_mask], db[labels[1]][final_mask], db[labels[2]][final_mask])
    ax.set_xlabel(labels[0])
    if limits[0] is not None:
        ax.set_xlim(*limits[0])
    ax.set_ylabel(labels[1])
    if limits[1] is not None:
        ax.set_ylim(*limits[1])
    ax.set_zlabel(labels[2])
    if limits[2] is not None:
        ax.set_zlim(*limits[2])

    annot = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords='offset points',
                bbox=dict(boxstyle='round', fc='w'),
                arrowprops=dict(arrowstyle='->'))
    annot.set_visible(False)
    #
    def update_annot(indices, xdata, ydata, zdata):
        x = numpy.mean(xdata[indices])
        y = numpy.mean(ydata[indices])
        z = numpy.mean(zdata[indices])
        names = db['Model'][indices]
        text = ', '.join(names)
        annot.set_text(text)
        x2, y2, _ = proj3d.proj_transform(x,y,z, ax.get_proj())
        annot.xy = x2, y2

    #
    def hover(event):
        if event.inaxes == ax:
            cont, ind = sp.contains(event)
            if cont:
                candidates = ind['ind']
                xdata = db[ax.get_xlabel()]
                ydata = db[ax.get_ylabel()]
                zdata = db[ax.get_zlabel()]
                update_annot(candidates, xdata, ydata, zdata)
                annot.set_visible(True)
            else:
                annot.set_visible(False)
            fig.canvas.draw()


    fig.canvas.mpl_connect('motion_notify_event', hover)
    plt.show()
    return ax, sp

def parse_parent_list():
    from collections import defaultdict
    model_to_parents=defaultdict(list)
    unique_parents = set()
    with open('parents.csv') as infile:
        lines = infile.read().split('\n')
    for line in lines[1:]:
        if not line:
            continue
        line = line.split(',')
        model_name = line[2]
        parent_id_and_chain = (line[3],line[4])
        model_to_parents[model_name].append(parent_id_and_chain)
        unique_parents.update([parent_id_and_chain])
    return model_to_parents, unique_parents

def find_parent(model, model_to_parents, templates):
    key = model.name[5:12]
    try:
        possible_templates = model_to_parents[key]
    except KeyError:
        return None
    for pt in possible_templates:
        template_name, template_chain = pt
        for t in templates:
            if t.name == template_name:
                if template_chain and t.chains[0].chain_id == template_chain:
                    return t
    return None


def load_possible_templates(session):
    import os
    model_to_parents, unique_parents = parse_parent_list()
    from chimerax.core.commands import open
    from chimerax.core.models import Model
    from chimerax.atomic import Residue
    import numpy
    templates = []
    # template_container = Model('templates', session)
    # session.models.add([template_container])
    # if not os.path.exists('templates'):
    #     os.mkdir('templates')
    for u_p in unique_parents:
        pdb_id, chain = u_p
        try:
            models = open.open(session, pdb_id)
            # If an ensemble, keep only the first model
            if len(models) > 1:
                session.models.close(models[1:])
            m = models[0]
            templates.append(m)
        except:
            continue
        m.residues[m.residues.polymer_types != Residue.PT_AMINO].delete()
        if chain != '':
            unique_chain_ids = numpy.unique(m.residues.chain_ids)
            if chain.lower() in unique_chain_ids and not chain.upper() in unique_chain_ids:
                chain = chain.lower()
            elif chain.upper() in unique_chain_ids and not chain.lower() in unique_chain_ids:
                chain = chain.upper()
            elif chain not in unique_chain_ids:
                # Throw up our hands. Probably a typo
                templates.remove(m)
                session.models.close([m])
                continue
                # raise RuntimeError("Chain {} not in PDB ID {}. Available chains are: {}".format(chain, pdb_id, ','.join(unique_chain_ids)))
            # print('Chain IDS: {}'.format(','.join(numpy.unique(m.residues.chain_ids))))
            # print('Deleting all residues with chain ID != {}'.format(chain))
            m.residues[m.residues.chain_ids != chain].delete()
        if len(numpy.unique(m.residues.chain_ids)) != 1:
            # The predictor hasn't bothered to specify the chain they actually
            # used as a template. This means we have more work to do.
            from chimerax.std_commands import split
            unaligned_chains = split.split_by_chain(m.atoms)
            unaligned_chains = [c[1] for c in unaligned_chains]
            # Break down into overlapping and keep the first instance of each.
            # Trust find_and_trim_matching_templates() to sort out which one is
            # the one we want
            keep_chains = []
            while len(unaligned_chains):
                ref_chain = unaligned_chains.pop(0)
                keep_chains.append(ref_chain)
                ref_chain = ref_chain[1]
                remaining_chains = []
                for c in unaligned_chains[1:]:
                    try:
                        alignment = align(session, ref_chain, c)
                    except:
                        remaining_chains.append(c)
                        continue
                    aligned_atoms = alignment[0]
                    if len(aligned_atoms) < 0.4*len(ref_chain.residues):
                        remaining_chains.append(c)
                unaligned_chains = remaining_chains
            # Remove the multi-chain model and add individual chain models
            new_models = [split.molecule_from_atoms(m, kc, name=m.name) for kc in keep_chains]
            session.models.add(new_models)
            templates.extend(new_models)
            templates.remove(m)
            session.models.close([m])

    # template_container.add(templates)
    return templates

def find_and_trim_matching_templates(session, target, templates, model_to_parents, cutoff_distance = 3,
        match_threshold = 0.4):
    from chimerax.core.commands import save as cxsave
    import os
    aligned_template_dirname = 'aligned_and_trimmed_templates'
    if not os.path.exists(aligned_template_dirname):
        os.mkdir(aligned_template_dirname)
    for t in reversed(templates):
        t.residues.ribbon_displays = False
        t.atoms[t.atoms.element_names == 'H'].delete()
        t.atoms.displays=True
        try:
            overall_result = align(session, target, t, cutoff_distance=None)
            best_result = align(session, target, t, cutoff_distance=cutoff_distance)
        except:
            templates.remove(t)
            try:
                print(t.id_string)
                session.models.close([t])
            except:
                print ('Failed on template {} before alignment'.format(t))
                return t
            continue
        aligned_atoms = best_result[0]
        if len(aligned_atoms) < len(target.residues)*match_threshold:
            templates.remove(t)
            try:
                print('Closing template {} due to too few residues overlapping with target'.format(t.id_string))
                session.models.close([t])
                continue
            except:
                print ('Failed on template {} after pruning'.format(t))
                return t
        # Delete everything in the template more than 20 residues outside the aligned region
        all_aligned_atoms = overall_result[0]
        aligned_residue_numbers = all_aligned_atoms.residues.numbers
        rmin = min(aligned_residue_numbers)-20
        rmax = max(aligned_residue_numbers)+20
        import numpy
        t_nums = t.residues.numbers
        delete_mask = numpy.logical_or(t_nums<rmin, t_nums>rmax)
        t.residues[delete_mask].delete()
        t.alignment_to_target = overall_result
        cxsave.save(session, os.path.join(aligned_template_dirname, t.name)+'.pdb', models=[t], rel_model=target)

logfile = open('align.log', 'wt')

def align(session, target, model, cutoff_distance=3, logfile=logfile):
    logfile.write('Aligning {} to {}\n'.format(model.name, target.name))
    logfile.flush()
    from chimerax.match_maker.match import match, defaults
    result = match(session, defaults['chain_pairing'], (target, [model]),
        defaults['matrix'], defaults['alignment_algorithm'],
        defaults['gap_open'], defaults['gap_extend'],
        cutoff_distance=cutoff_distance, always_raise_errors=True)[0]

    # result = cmd_match(session, model.atoms, to=target.atoms,
    # cutoff_distance=cutoff_distance, always_raise_errors=True)[0]

    # Returned arrays of aligned atoms are paired, but are in random order with
    # respect to the actual models. Re-sort them before returning
    import numpy
    model_atoms = result[0]
    target_atoms = result[1]
    sort_order = numpy.argsort(target_atoms.residues.numbers)
    model_atoms = model_atoms[sort_order]
    target_atoms = target_atoms[sort_order]
    return (model_atoms, target_atoms, *result[2:])


def sidechain_buried_score(residue):
    '''
    Defines how "buried" a sidechain is by counting the number of heavy atoms
    from other residues coming within 4A of any heavy atom from the sidechain.
    The returned score is a value ranging from 0 to 1, where 0 indicates no
    contact with other atoms, and 1 indicates 3 or more other atoms per
    sidechain atoms. The score scales linearly with the number of contacting
    atoms in between these values.
    '''
    from chimerax.core.geometry import find_close_points
    from chimerax.atomic import Residues
    import numpy
    r = residue
    m = r.structure
    other_residues = m.residues.subtract(Residues([r]))
    sidechain_atoms = r.atoms[numpy.logical_not(numpy.in1d(r.atoms.names, ['N', 'C', 'CA', 'O']))]
    if not len(sidechain_atoms):
        return 0
    other_atoms = other_residues.atoms
    cp = find_close_points(sidechain_atoms.coords, other_atoms.coords, 4.0)[1]
    score = (len(cp)/len(sidechain_atoms))/3
    if score > 1:
        score = 1
    return score


def chord_formula(angle_deltas):
    from numpy import cos
    return (1-cos(angle_deltas))/2

def backbone_score(phi_deltas, psi_deltas, chi_deltas):
    import numpy
    combined = numpy.column_stack([chord_formula(d) for d in (phi_deltas, psi_deltas, chi_deltas)])
    return numpy.nanmean(combined, axis=1)

def sidechain_score(burial_scores, n_chis, chi1_deltas, chi2_deltas, chi1_cutoff = radians(30)):
    import numpy
    no_chi1 = numpy.isnan(chi1_deltas)
    no_chi2 = numpy.isnan(chi2_deltas)
    scores = numpy.zeros(chi1_deltas.shape)
    # Sidechains with no rotamers do not contribute to the score
    scores[n_chis==0] = numpy.nan
    # If a rotameric sidechain is present in the target but not in the model, it
    # gets the maximum possible score
    scores[numpy.logical_and(n_chis>0, no_chi1)] = 1
    chi1_scores = (1+numpy.cos(chi1_deltas))/2
    chi2_scores = (1+numpy.cos(chi2_deltas))/2
    single_chi_filter = numpy.logical_and(~no_chi1, no_chi2)
    two_chi_filter = numpy.logical_and(~no_chi1, ~no_chi2)

    # print('Scores: {}\nBurial_scores: {}\nChi1_scores: {}'.format(scores, burial_scores, chi1_scores))
    scores[single_chi_filter] = (burial_scores[single_chi_filter] *
        (1-chi1_scores[single_chi_filter]))

    scores[two_chi_filter] = (burial_scores[two_chi_filter] *
        (1-2/3*chi1_scores[two_chi_filter]
          - 1/3 * numpy.exp(-(chi1_deltas[two_chi_filter]/chi1_cutoff)**2)
                * chi2_scores[two_chi_filter]))
    return scores


def compare_torsions(session, target, decoy):
    '''
    Residue-by-residue comparison of key torsion angles for two models.
    target and decoy should be two models of the same protein, each containing
    a single chain.
    '''
    import numpy
    from chimerax.match_maker import match
    from chimerax.match_maker.settings import defaults
    from chimerax.atomic import Residues
    from chimerax.isolde import session_extensions as sx
    from math import pi, sqrt
    proper_dihedral_mgr = sx.get_proper_dihedral_mgr(session)
    protein_dihedral_dict = proper_dihedral_mgr.dihedral_dict['residues']['protein']

    backbone_dihedral_names = ('phi', 'psi', 'omega')
    sidechain_dihedral_names = ('chi1', 'chi2', 'chi3', 'chi4')
    t_chain_ids = target.chains.chain_ids
    d_chain_ids = decoy.chains.chain_ids
    if len(t_chain_ids) !=1 or len(d_chain_ids) !=1:
        raise RuntimeError('Target and decoy must each have a single chain!')
    if hasattr(decoy, 'alignment_to_target'):
        # Alignment has already been determined. re-use it
        alignment_result = decoy.alignment_to_target
    else:
        # Align once with no cut-off for maximum overlap
        try:
            alignment_result = decoy.alignment_to_target = align(session, target, decoy, cutoff_distance=None)
        except:
            raise
            print('Failed to align {} to {} - bailing out on this analysis!'.format(decoy, target))
            return None
        # Align again with default distance cutoff for best alignment
        align(session, target, decoy)

    d_ca = alignment_result[0]
    t_ca = alignment_result[1]

    ca_ca_distances = numpy.linalg.norm(d_ca.scene_coords-t_ca.scene_coords, axis=1)


    d_res = d_ca.residues
    t_res = t_ca.residues

    residue_numbers = t_res.numbers
    backbone_mean_error = numpy.empty(len(t_res))
    sidechain_mean_error = numpy.empty(len(t_res))
    twists = numpy.zeros(len(t_res), numpy.bool)
    cis_trans_flips = numpy.zeros(len(t_res), numpy.bool)
    target_num_chi = numpy.zeros(len(t_res), int)
    deltas = {
    'phi': numpy.ones(len(t_res))*numpy.nan,
    'psi': numpy.ones(len(t_res))*numpy.nan,
    'omega': numpy.ones(len(t_res))*numpy.nan,
    'chi1': numpy.ones(len(t_res))*numpy.nan,
    'chi2': numpy.ones(len(t_res))*numpy.nan,
    'chi3': numpy.ones(len(t_res))*numpy.nan,
    'chi4': numpy.ones(len(t_res))*numpy.nan
    }
    sidechain_burial_scores = numpy.ones(len(t_res))*numpy.nan

    for i, (tr, dr) in enumerate(zip(t_res, d_res)):
        nchi = protein_dihedral_dict[tr.name]['nchi']
        target_num_chi[i] = nchi
        symm = protein_dihedral_dict[tr.name]['symm']
        t_name = tr.name
        d_name = dr.name
        chi_adjustment = None
        if t_name != d_name:
            if t_name in _torsion_adjustments.keys():
                chi_adjustment = _torsion_adjustments[t_name]
            elif d_name in _torsion_adjustments.keys():
                chi_adjustment = list(_torsion_adjustments[d_name])
                chi_adjustment[1] = -chi_adjustment[1]

        abs_b_mean = 0
        for bd_name in backbone_dihedral_names:
            td = proper_dihedral_mgr.get_dihedral(tr, bd_name)
            dd = proper_dihedral_mgr.get_dihedral(dr, bd_name)
            if bd_name == 'omega' and dd is not None:
                dd_angle = dd.angle
                if abs(dd_angle) > TWISTED_THRESHOLD and abs(dd_angle) < pi-TWISTED_THRESHOLD:
                    twists[i] = True
                if td is not None:
                    if abs(minimum_difference(dd_angle, td.angle)) > pi/2:
                        cis_trans_flips[i] = True
            if td is None or dd is None:
                # No penalty for missing backbone dihedrals (only occur at
                # chain breaks)
                abs_b_mean += 0
            else:
                diff = minimum_difference(td.angle, dd.angle)
                deltas[bd_name][i] = diff
                abs_b_mean += abs(diff)
        abs_b_mean = abs_b_mean/3
        backbone_mean_error[i] = abs_b_mean

        abs_s_mean = 0
        sum_weights = 0

        # If a sidechain has no rotamers, it doesn't contribute to the score
        buried_score = sidechain_buried_score(tr)
        sidechain_burial_scores[i] = buried_score

        if nchi==0:
            sidechain_mean_error[i] = numpy.nan
            continue

        if buried_score == 0:
            sidechain_mean_error[i] = 0
            continue

        for j in range(nchi):
            sd_name =sidechain_dihedral_names[j]
            td = proper_dihedral_mgr.get_dihedral(tr, sd_name)
            dd = proper_dihedral_mgr.get_dihedral(dr, sd_name)

            if td is not None:
                # downweight angles further from the backbone
                weight = 1/2**j
                sum_weights += weight
                if dd is None:
                    # Big penalty for missing sidechain in decoy
                    abs_s_mean += pi * weight
                else:
                    td_angle = td.angle
                    if chi_adjustment is not None:
                        if sd_name == chi_adjustment[0]:
                            td_angle += chi_adjustment[1]
                    dd_angle = dd.angle
                    if symm and j == nchi-1:
                        if td_angle < 0:
                            td_angle += pi
                        if dd_angle < 0:
                            dd_angle += pi
                    diff = minimum_difference(td_angle, dd_angle)
                    deltas[sd_name][i] = diff
                    abs_s_mean += abs(diff) * weight
        if sum_weights > 0:
            abs_s_mean = abs_s_mean/sum_weights
            sidechain_mean_error[i]=abs_s_mean*buried_score
        else:
            sidechain_mean_error[i] = numpy.nan
    backbone_mean_error = numpy.degrees(backbone_mean_error)
    sidechain_mean_error = numpy.degrees(sidechain_mean_error)
    existing_sidechain_mask = numpy.logical_not(numpy.isnan(sidechain_mean_error))
    results = {
        'Target residue numbers': residue_numbers,
        'CA-CA distances': ca_ca_distances,
        'Backbone angle error': backbone_mean_error,
        'Overall backbone angle error': numpy.mean(backbone_mean_error),
        'Weighted chi angle error': sidechain_mean_error,
        'Overall weighted chi angle error': numpy.mean(sidechain_mean_error[existing_sidechain_mask]),
        'Twisted peptide bonds': twists,
        'Peptide cis/trans flips': cis_trans_flips,
        'Aligned residues': (t_res, d_res),
        'Torsion deltas': deltas,
        'Backbone score': backbone_score(deltas['phi'], deltas['psi'], deltas['omega']),
        'Sidechain score': sidechain_score(sidechain_burial_scores, target_num_chi, deltas['chi1'], deltas['chi2'])
        }
    return results


def model_vs_target(session, target, model, template=None, plot=True):
    models = (target, template, model)
    for i, m in enumerate(models):
        if m is not None and len(m.chains.chain_ids) != 1:
            if i < 2:
                if template is not None:
                    raise TypeError("i={}; target={}:{}, model={}:{}, template={}:{}".format(i, target.name, target.id_string, model.name, model.id_string, template.name, template.id_string))
                else:
                    raise TypeError("i={}; target={}:{}, model={}:{}, template=None".format(i, target.name, target.id_string, model.name, model.id_string))

                raise TypeError("Model {}:{} contains {} chains!".format(m.name, m.id_string, len(m.chains)))
            else:
                # Something weird is wrong with the template. Just throw it out.
                template = None
    if template is not None:
        template_vs_target = compare_torsions(session, target, template)
        template_name = template.name
    else:
        template_vs_target = None
        template_name = "N/A"
    print('Model: {} Target: {}'.format(model.name, target.name))
    model_vs_target = compare_torsions(session, target, model)
    if model_vs_target is None:
        return None

    if plot:

        fig = plot_comparison("Model {} torsion-angle analysis (template = {})".format(os.path.splitext(model.name)[0], template_name),
            "Model vs Target", model_vs_target,
            "Template vs Target", template_vs_target
            )
    else:
        fig = None
    return (fig, template_vs_target, model_vs_target)


def minimum_difference(angle1, angle2):
    from math import pi
    a = angle1-angle2
    return abs(a+pi)%(2*pi)-pi




def plot_results(results_dict, key):
    import numpy
    from matplotlib import pyplot as plt
    residue_numbers = results_dict[key]['Target residue numbers']
    backbone_rms = results_dict[key]['Backbone angle error']
    sidechain_rms = results_dict[key]['Weighted chi angle error']

    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(residue_numbers, backbone_rms)
    ax1.set_title("Backbone dihedral mean error")
    ax1.set_ylabel("mean abs angle difference")

    ax2.plot(residue_numbers, sidechain_rms)
    ax2.set_title("Weighted sidechain dihedral error")
    ax2.set_ylabel("mean abs angle difference")

    plt.show()

def process_directories(session, paths, reprocess=False):
    import os
    from glob import glob
    for p in paths:
        pflag = glob(os.path.join(p, 'processed'))
        if len(pflag):
            if not reprocess:
                continue
            os.remove(pflag[0])
        print('Processing {}'.format(os.path.basename(p)), flush=True)
        process_directory(session, p)
        session.models.close(session.models.list())

def process_directory(session, full_path):
    print('Processing {}'.format(full_path))
    import os
    from matplotlib import pyplot as plt
    import gc
    target_name = os.path.basename(full_path)
    os.chdir(full_path)
    # Close all currently open models
    session.models.close(session.models.list())
    try:
        target = load_target(session, target_name)
    except:
        print("Target {} missing!".format(target_name))
        return
    model_to_parents, unique_parents = parse_parent_list()
    templates = load_possible_templates(session)
    find_and_trim_matching_templates(session, target, templates, model_to_parents)
    with open(target_name+'_torsion_comparison.csv', 'wt') as outfile:
        write_outfile_header(outfile)
        from glob import glob
        from chimerax.core.commands import open as cxopen
        all_model_files = glob('*.pdb')
        for mfile in all_model_files:
            m = cxopen.open(session, mfile)[0]
            parent = find_parent(m, model_to_parents, templates)
            result = model_vs_target(session, target, m, template=parent)
            if result is None:
                session.models.close([m])
                continue

            fig, t_vs_t, m_vs_t = result

            delta_dict = {}
            delta_dict['Model'] = {key: vals.tolist() for key, vals in m_vs_t['Torsion deltas'].items()}
            if parent is not None:
                delta_dict['Template'] = {key: vals.tolist() for key, vals in t_vs_t['Torsion deltas'].items()}
            import json
            with open(os.path.join(full_path, mfile+'torsion_deltas.json'), 'wt') as jsonfile:
                json.dump(delta_dict, jsonfile)

            model_name = os.path.splitext(m.name)[0]
            if parent is None:
                parent_name = "none"
            else:
                parent_name = parent.name
            fig.savefig(model_name+'_vs_'+parent_name+'_torsion_comparison_new.png')
            fig.clf()
            plt.close(fig)
            gc.collect()
            append_comparison_to_file(outfile, model_name, m_vs_t, parent_name, t_vs_t)
            session.models.close([m])
    with open('processed', 'wt') as pflag:
        pflag.write('1')

def process_refinement_directories(session, paths, reprocess=False):
    import os
    from glob import glob
    for p in paths:
        pflag = glob(os.path.join(p, 'processed'))
        if len(pflag):
            if not reprocess:
                continue
            os.remove(pflag[0])
        print('Processing {}'.format(os.path.basename(p)), flush=True)
        process_refinement_directory(session, p)

def process_refinement_directory(session, full_path):
    import os
    import gc
    from matplotlib import pyplot as plt
    session.models.close(session.models.list())
    refinement_name = os.path.basename(full_path)
    os.chdir(full_path)
    target = load_target(session, refinement_name)
    starting_model_dirname, starting_model_filename = refinement_starting_models[refinement_name]
    from chimerax.core.commands import open as cxopen
    starting_model = cxopen.open(session, os.path.join(_prediction_dir, starting_model_dirname, starting_model_filename))[0]
    starting_model_name = os.path.splitext(starting_model.name)[0]
    with open(refinement_name+'_torsion_comparison.csv', 'wt') as outfile:
        write_outfile_header(outfile)
        # First do the starting model
        result = model_vs_target(session, target, starting_model, plot=False)
        _, t_vs_t, m_vs_t = result
        append_comparison_to_file(outfile, 'starting_model', m_vs_t, 'none', t_vs_t)
        from glob import glob
        from chimerax.core.commands import open as cxopen
        all_model_files = glob('*.pdb')
        for mfile in all_model_files:
            m = cxopen.open(session, mfile)[0]
            result = model_vs_target(session, target, m, template=starting_model)
            if result is None:
                session.models.close([m])
                continue
            fig, t_vs_t, m_vs_t = result
            model_name = os.path.splitext(m.name)[0]
            fig.savefig(model_name+'_vs_'+starting_model_name+'_torsion_comparison.png')
            plt.close(fig)
            gc.collect()
            append_comparison_to_file(outfile, model_name, m_vs_t, starting_model_name, t_vs_t)
            session.models.close([m])
    with open('processed', 'wt') as pflag:
        pflag.write('1')




def write_outfile_header(output_file):
    base_entries = ['Name', 'Coverage', 'Average Backbone Error', 'Backbone Score',
        'Twisted Peptide Count', 'Twisted Peptide Frac', 'cis/trans Flips',
        'cis/trans Flip Frac', 'Average Sidechain Error', 'Sidechain Score']

    header_entries = (['Model ' + e for e in base_entries] +
                     ['Template ' + e for e in base_entries])
    # header_entries = ['Model Name', 'Model Coverage', 'Model Average Backbone Error',
    #     'Model Twisted Peptide Count', 'Model Twisted Peptide Frac',
    #     'Model cis/trans Flips', 'Model cis/trans Flip Frac',
    #     'Model Average Sidechain Error', 'Template Name', 'Template Coverage',
    #     'Template Average Backbone Error', 'Template Twisted Peptide Count',
    #     'Template Twisted Peptide Frac', 'Template cis/trans Flips',
    #     'Template cis/trans Flip Frac', 'Template Average Sidechain Error']
    output_file.write(','.join(header_entries)+'\n')

def append_comparison_to_file(outfile, model_name, model_result, template_name, template_result):
    '''
    Append results to a comma-separated text file
    '''
    import numpy
    line_entry = []
    line_entry.append(model_name) # Model Name
    target_res, model_res = model_result['Aligned residues']
    target = target_res.unique_structures[0]
    model = model_res.unique_structures[0]
    line_entry.append('{:0.3f}'.format(len(target_res)/len(target.residues))) # Model Coverage
    line_entry.append('{:0.3f}'.format(model_result['Overall backbone angle error'])) # Model Average Backbone Error
    line_entry.append('{:0.3f}'.format(numpy.nanmean(model_result['Backbone score']))) # Chord-based score on phi/psi/omega
    model_twisted_peptide_count = sum(model_result['Twisted peptide bonds'])
    line_entry.append(str(model_twisted_peptide_count)) # Model Twisted Peptide Count
    line_entry.append('{:.3f}'.format(model_twisted_peptide_count/len(model.residues))) # Frac
    model_cis_trans_flips = sum(model_result['Peptide cis/trans flips'])
    line_entry.append(str(model_cis_trans_flips)) # Model cis/trans Flips
    line_entry.append('{:.3f}'.format(model_cis_trans_flips/len(model.residues))) # Frac
    line_entry.append('{:.3f}'.format(model_result['Overall weighted chi angle error'])) # Model Average Sidechain Error
    line_entry.append('{:.3f}'.format(numpy.nanmean(model_result['Sidechain score']))) # Chord-based score on chi1/chi2
    if template_result is not None:
        line_entry.append(template_name) # Template Name
        target_res, template_res = template_result['Aligned residues']
        template = template_res.unique_structures[0]
        line_entry.append('{:0.3f}'.format(len(target_res)/len(target.residues))) # Template Coverage
        line_entry.append('{:0.3f}'.format(template_result['Overall backbone angle error'])) # Template Average Backbone Error
        line_entry.append('{:0.3f}'.format(numpy.nanmean(template_result['Backbone score']))) # Chord-based score on phi/psi/omega
        template_twisted_peptide_count = sum(template_result['Twisted peptide bonds'])
        line_entry.append(str(template_twisted_peptide_count)) # Template Twisted Peptide Count
        line_entry.append('{:.3f}'.format(template_twisted_peptide_count/len(template.residues))) # Frac
        template_cis_trans_flips = sum(template_result['Peptide cis/trans flips'])
        line_entry.append(str(template_cis_trans_flips)) # Template cis/trans Flips
        line_entry.append('{:.3f}'.format(template_cis_trans_flips/len(template.residues))) # Frac
        line_entry.append('{:.3f}'.format(template_result['Overall weighted chi angle error'])) # Template Average Sidechain Error
        line_entry.append('{:.3f}'.format(numpy.nanmean(template_result['Sidechain score']))) # Chord-based score on chi1/chi2
    else:
        line_entry.extend(['N/A']*10)
    outfile.write(','.join(line_entry)+'\n')


def plot_comparison(name, model_name, model_result, template_name, template_result, fig=None):
    from matplotlib import pyplot as plt
    from matplotlib import lines as mlines
    import numpy

    fig, (ax1, ax1_delta, ax2, ax2_delta) = plt.subplots(4,1, figsize=(10,8))
    fig.suptitle(name)
    res1 = model_result['Target residue numbers']
    br1 = model_result['Backbone score']
    sr1 = model_result['Sidechain score']
    # br1 = model_result['Backbone angle error']
    # sr1 = model_result['Weighted chi angle error']
    tw1 = model_result['Twisted peptide bonds']
    ct1 = model_result['Peptide cis/trans flips']
    m_ca_distances = model_result['CA-CA distances']
    m_bb_mean = numpy.nanmean(br1)
    m_sc_mean = numpy.nanmean(sr1)
    # m_bb_mean = model_result['Overall backbone angle error']
    # m_sc_mean = model_result['Overall weighted chi angle error']

    if template_result is not None:
        res2 = template_result['Target residue numbers']
        br2 = template_result['Backbone score']
        sr2 = template_result['Sidechain score']
        # br2 = template_result['Backbone angle error']
        # sr2 = template_result['Weighted chi angle error']
        tw2 = template_result['Twisted peptide bonds']
        ct2 = template_result['Peptide cis/trans flips']
        t_ca_distances = template_result['CA-CA distances']
        t_bb_mean = numpy.nanmean(br2)
        t_sc_mean = numpy.nanmean(sr2)
        # t_bb_mean = template_result['Overall backbone angle error']
        # t_sc_mean = template_result['Overall weighted chi angle error']

    else:
        res2 = br2 = sr2 = tw2 = ct2 = t_bb_mean = t_sc_mean = None

    import numpy


    res1_filled = numpy.arange(res1.min(), res1.max()+1)
    mask = numpy.in1d(res1_filled, res1)

    ca1_filled_vals = numpy.ones(len(res1_filled))*numpy.nan
    ca1_filled_points = numpy.ones(len(res1_filled))*numpy.nan
    ca1_filled_points[mask] = -0.5
    ca1_filled_vals[mask] = m_ca_distances

    br1_filled = numpy.ones(len(res1_filled))*numpy.nan
    br1_filled[mask] = br1
    sr1_filled = numpy.ones(len(res1_filled))*numpy.nan
    sr1_filled[mask] = sr1
    tw1_filled = numpy.ones(len(res1_filled))*numpy.nan
    tw1_filled[mask] = tw1
    tw1_filled[tw1_filled==0] = numpy.nan
    ct1_filled = numpy.ones(len(res1_filled))*numpy.nan
    ct1_filled[mask] = ct1

    ct1_filled[ct1_filled==0] = numpy.nan

    if template_result is not None:
        res2_filled = numpy.arange(res2.min(), res2.max()+1)
        mask = numpy.in1d(res2_filled, res2)

        ca2_filled_vals = numpy.ones(len(res2_filled))*numpy.nan
        ca2_filled_points = numpy.ones(len(res2_filled))*numpy.nan
        ca2_filled_points[mask] = 0.5
        ca2_filled_vals[mask] = t_ca_distances


        br2_filled = numpy.ones(len(res2_filled))*numpy.nan
        br2_filled[mask] = br2
        tw2_filled = numpy.ones(len(res2_filled))*numpy.nan
        tw2_filled[mask] = tw2
        tw2_filled[tw2_filled==0] = numpy.nan
        ct2_filled = numpy.ones(len(res2_filled))*numpy.nan
        ct2_filled[mask] = ct2
        ct2_filled[ct2_filled==0] = numpy.nan
        sr2_filled = numpy.ones(len(res2_filled))*numpy.nan
        sr2_filled[mask] = sr2

        res_combined = numpy.intersect1d(res1, res2)
        res_filled = numpy.arange(res_combined.min(), res_combined.max()+1)
        res_combined_filled = numpy.arange(res_combined.min(), res_combined.max()+1)
        mask_combined_filled = numpy.in1d(res_filled, res_combined)

        mask_filled_1 = numpy.in1d(res1, res_combined)
        br1_combined_filled = numpy.ones(len(res_combined_filled))*numpy.nan
        br1_combined_filled[mask_combined_filled] = br1[mask_filled_1]
        sr1_combined_filled = numpy.ones(len(res_combined_filled))*numpy.nan
        sr1_combined_filled[mask_combined_filled] = sr1[mask_filled_1]

        mask_filled_2 = numpy.in1d(res2, res_combined)
        br2_combined_filled = numpy.ones(len(res_combined_filled))*numpy.nan
        br2_combined_filled[mask_combined_filled] = br2[mask_filled_2]
        sr2_combined_filled = numpy.ones(len(res_combined_filled))*numpy.nan
        sr2_combined_filled[mask_combined_filled] = sr2[mask_filled_2]


    if template_result is None:
        template_backbone_mean = 'N/A'
        template_twisted_frac = 'N/A'
        template_cis_trans = 'N/A'
        template_sidechain_mean = 'N/A'
    else:
        template_backbone_mean = '{:.2f}'.format(t_bb_mean)
        template_twisted_frac = '{:.4f}'.format(sum(tw2)/len(tw2))
        template_cis_trans = '{}'.format(sum(ct2))
        template_sidechain_mean = '{:.2f}'.format(t_sc_mean)

    lines = []
    labels = []

    if template_result is None:
        abs_x_min = min(res1_filled)
        abs_x_max = max(res1_filled)
    else:
        abs_x_min = min((min(res1_filled), min(res2_filled)))
        abs_x_max = max((max(res1_filled), max(res2_filled)))
    abs_x_range = abs_x_max-abs_x_min
    padding = 0.025*abs_x_range
    x_axis_limits = [abs_x_min-padding, abs_x_max+padding]

    ax1.set_title("Backbone torsion penalty (template mean={}, model mean={:.2f})".format(template_backbone_mean, m_bb_mean))
    # ax1.set_ylabel("(|𝚫ω|+|𝚫φ|+|𝚫ψ|) / 3")
    ax1.set_ylabel("Backbone penalty")
    #ax1.set_ylabel("mean abs angle difference")
    if template_result is not None:
        ax1.plot(res2_filled, br2_filled, 'blue')
        lines.append(mlines.Line2D([],[], color='blue'))
        labels.append("Template")
    ax1.plot(res1_filled, br1_filled, 'red')

    lines.append(mlines.Line2D([],[],color='red'))
    labels.append("Model")
    ax1.set_ylim([0, 1])
    ax1.set_xlim(x_axis_limits)
    ax1.text(min(res1_filled), 0.889, "Template twisted peptide bond fraction: {}".format(template_twisted_frac))
    ax1.text(min(res1_filled), 0.722, "Model twisted peptide bond fraction: {:0.4f}".format(sum(tw1)/len(tw1)))
    ax1.text(max(res1_filled), 0.889, "Template cis<->trans flips: {}".format(template_cis_trans), horizontalalignment='right')
    ax1.text(max(res1_filled), 0.722, "Model cis<->trans flips: {}".format(sum(ct1)), horizontalalignment='right')



    from matplotlib.collections import LineCollection
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    colorbar_label_x = 0.89*(x_axis_limits[1]-x_axis_limits[0]) + x_axis_limits[0]

    ax1_delta.set_title("Model score - Template score")
    if template_result is not None:
        t_pts = numpy.array([res2_filled, ca2_filled_points]).T.reshape(-1,1,2)
        t_segs = numpy.concatenate([t_pts[:-1], t_pts[1:]], axis=1)
        norm = plt.Normalize(0, 5)
        t_lc = LineCollection(t_segs, cmap='Greens', norm=norm)
        t_lc.set_array(ca2_filled_vals)
        t_lc.set_linewidth(45)
        t_lc.set_alpha(0.3)
        t_line = ax1_delta.add_collection(t_lc)

        ax1_delta.text(colorbar_label_x, 0.7222, "Template CA offset (Å)", horizontalalignment='right')
        t_cbaxes = inset_axes(ax1_delta, width="10%", height="7.5%", loc='upper right')
        plt.colorbar(mappable=t_line, cax=t_cbaxes, orientation="horizontal", ticks=[0., 5])


        ax1_delta.plot(res_combined_filled, br1_combined_filled-br2_combined_filled, 'black')
        ax1_delta.plot(res2_filled, tw2_filled.astype(float)*0.667, 'b^') # twisted in template
        ax1_delta.plot(res2_filled, ct2_filled.astype(float)*0.556, 'bx') # cis<->trans between template and target

    m_pts = numpy.array([res1_filled, ca1_filled_points]).T.reshape(-1,1,2)
    m_segs = numpy.concatenate([m_pts[:-1], m_pts[1:]], axis=1)
    norm = plt.Normalize(0, 5)
    m_lc = LineCollection(m_segs, cmap='Purples', norm=norm)
    m_lc.set_array(ca1_filled_vals)
    m_lc.set_linewidth(45)
    m_lc.set_alpha(0.3)
    m_line = ax1_delta.add_collection(m_lc)

    ax1_delta.set_ylabel('Delta')
    ax1_delta.text(colorbar_label_x, -0.9167, "Model CA offset (Å)", horizontalalignment='right')
    m_cbaxes = inset_axes(ax1_delta, width="10%", height="7.5%", loc='lower right')
    plt.colorbar(mappable=m_line, cax=m_cbaxes, orientation="horizontal", ticks=[0., 5])

    #fig.colorbar()



    ax1_delta.plot(res1_filled, tw1_filled.astype(float)*-0.667, 'r^') # twisted in model
    lines.append(mlines.Line2D([],[],color='black', marker='^', linestyle='None'))
    labels.append("ω twisted >30° from plane")
    ax1_delta.plot(res1_filled, ct1_filled.astype(float)*-0.556, 'rx') # cis<->trans between model and target
    lines.append(mlines.Line2D([],[],color='black', marker='x', linestyle='None'))
    labels.append("ω cis ↔ trans")
    ax1_delta.set_ylim([-1,1])
    ax1_delta.set_xlim(x_axis_limits)

    ax2.set_title("Sidechain torsion penalty (template mean={}, model mean={:.2f})".format(template_sidechain_mean, m_sc_mean))
    # ax2.set_ylabel("β𝚺(w|𝚫𝜒|) / 𝚺w")
    ax2.set_ylabel('Sidechain penalty')
    # ax2.set_ylabel("mean abs angle difference")
    if template_result is not None:
        ax2.plot(res2_filled, sr2_filled, 'blue')
    ax2.plot(res1_filled, sr1_filled, 'red')
    ax2.set_ylim([0,1])
    ax2.set_xlim(x_axis_limits)

    ax2_delta.set_title("Model score - Template score")
    ax2_delta.set_ylabel("Delta")
    if template_result is not None:
        ax2_delta.plot(res_combined_filled, sr1_combined_filled-sr2_combined_filled, 'black')
        ax2_delta.set_xlim([res1_filled.min(), res1_filled.max()])
    ax2_delta.set_ylim([-1,1])
    ax2_delta.set_xlim(x_axis_limits)

    fig.legend(lines, labels, 'upper right')
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    plt.show()
    return fig


def _char_in_font(unicode_char, font):
    for cmap in font['cmap'].tables:
        if cmap.isUnicode():
            if ord(unicode_char) in cmap.cmap:
                return True
    return False

def _find_suitable_fonts(chars):
    '''
    Returns a list of fonts containing all the necessary Unicode symbols.
    '''
    from fontTools.ttLib import TTFont
    import matplotlib.font_manager as mfm
    fonts = []
    font_info = [(f.fname, f.name) for f in mfm.fontManager.ttflist]
    for i, font in enumerate(font_info):
        missing = False
        f = TTFont(font[0])
        for c in chars:
            if not _char_in_font(c, f):
                missing = True
                break
        if not missing:
            fonts.append(font)
    return fonts

def process_easy_targets(session, reprocess=False, start=0, end=None):
    import os
    if end is None or end > len(easy_targets):
        end = len(easy_targets)
    print('Easy targets: {}'.format(easy_targets))
    process_directories(session, [os.path.join(_prediction_dir, d) for d in easy_targets[start:end]], reprocess=reprocess)

def process_hard_targets(session, reprocess=False):
    import os
    process_directories(session, [os.path.join(_prediction_dir, d) for d in hard_targets], reprocess=reprocess)
