from collections import defaultdict
from itertools import chain
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

import numpy as np
from toposort import toposort
from mzqc import MZQCFile as mzqc
import pyopenms as oms

def getLongestTag(spec: oms.MSSpectrum, aa_weights: List[float], tol: float=0.5) -> int:
    # TODO spec.getPrecursors()[0].getCharge() > 2 consider aa_weights for doubly charged or always?
    # what about internal fragments and modifications
    # aa_weights_z1 = np.array(list({r.getMonoWeight(1) for r in oms.ResidueDB().getResidues('AllNatural')}),dtype=float)
    # aa_weights_z2 = np.array(list({r.getMonoWeight(2) for r in oms.ResidueDB().getResidues('AllNatural')}),dtype=float)/2
    if spec.getMSLevel() == 1:
        return -1
    if not spec.isSorted():
        spec.sortByPosition()

    edges: List[Tuple[int,int]] = list() 
    node_dependencies: Dict[int,Set[Any]] = defaultdict(set)
    path_score: Dict[int,int] = defaultdict(lambda: np.NINF)
    for i in range(0,spec.size()-1):
        for j in range(i,spec.size()):
            dist = spec[j].getMZ()-spec[i].getMZ()
            if np.any(np.isclose(dist,aa_weights, atol=tol)):
                edges.append((i,j))
                node_dependencies[j].add(i)
    topological = list(toposort(node_dependencies))
    if len(topological)<=1:
        return 0
    for obound in topological[0]:
        path_score[obound] = 0
    # topological = list(toposort({2: {11},
    #         9: {11, 8, 10},
    #         10: {11, 3},
    #         11: {7, 5},
    #         8: {7, 3},
    # }))
    # edges = [(3,8),(3,10),(5,11),(7,8),(7,11),(8,9),(11,2),(11,9),(11,10),(10,9)]
    edge_array = np.array(edges, ndmin = 2)
    edge_sort_order = list(chain.from_iterable(topological))
    for u in edge_sort_order:
        for edge in edge_array[edge_array[:,0] == u]:
            if path_score[edge[1]] < path_score[edge[0]] + 1:  # edgecost always 1
                path_score[edge[1]] = path_score[edge[0]] + 1

    return max(path_score.values())

def getMassTraceMatchingMS2(exp: oms.MSExperiment, tol: float=0.5) -> List[mzqc.QualityMetric]:
    mts: List[oms.MassTrace] = list()
    oms.MassTraceDetection().run(exp,mts,0)  # since 2.5.0 with 3rd argument
    mts_coord = np.array([[m.getCentroidMZ(),m.getCentroidRT()] for m in mts])
    # ms2_coord = np.array([[s.getPrecursors()[0].getMZ(), s.getRT()] for s in exp if s.getMSLevel()==2])
    
    for s in exp:
        if s.getMSLevel()==2:
            mz_matches = np.isclose(mts_coord[:,0], s.getPrecursors()[0].getMZ(), atol=tol)
            rt_dist_per_match = np.abs(mts_coord[np.where(mz_matches)][:,1] - s.getRT())
            match_idx_in_dist = np.argwhere(mz_matches)  # indices of match only in mts and mts_coord 
            closest_rt_rowidx = rt_dist_per_match.argmin()  # index in match_only distances array
            # rt_dist_per_match[closest_rt_rowidx] == mts_coord[match_idx[closest_rt_rowidx][0]][1]-s.getRT()
            closest_match_mt = mts[match_idx_in_dist[closest_rt_rowidx][0]]

            np.partition(rt_dist_per_match,2)[2-1]  # 2nd closest dist
            np.partition(rt_dist_per_match,1)[1-1]  # closest dist
            closest_match_mt.getSize()  # peaks
            closest_match_mt.getTraceLength()  # seconds
            closest_match_mt.getFWHM()  # seconds - what if 0 or getTraceLength()?
            closest_match_mt.getMaxIntensity(False)

            # NB precursor intensity is always 0! 
            # NB masstrace does not store peak intensities (except max and sum)
            # 4 categories for MS2 regarding sampling 
            # -2 (out of trace, before centr RT) ; -1 (in trace, before centr RT) ;1 (in trace, after centr RT) ;2 (out of trace, after centr RT) ;
            rt_1st = np.min(closest_match_mt.getConvexhull().getHullPoints()[:,0])
            rt_last = np.max(closest_match_mt.getConvexhull().getHullPoints()[:,0])
            rt_centr = closest_match_mt.getCentroidRT()
            # np.digitize(s.getRT(),[rt_1st,rt_centr,rt_last])
            if s.getRT() > rt_centr:  # 'after' categ
                if s.getRT() > rt_last:
                    return 2
                else:
                    return 1
            else:  # 'before' categ
                if s.getRT() < rt_1st:
                    return -2
                else:
                    return -1


    # get mts
    # for each ms2 find mz matched within tol
    #     pick closest match in RT np.array([1.1,2.2,3.3,.8,5.5,6.6,.7,8.8]
    #     report:
    #         how close is the closest
    #         on which side is the precursor
    #         how wobbly is the mt (mz sd)
    #         how long is the mt (fwhm)
    #         are there others close? (next closest)	


    # computePeakArea(
    # computeSmoothedPeakArea(
    # findMaxByIntPeak(
    # estimateFWHM(
    # getFWHM(
    # getSmoothedIntensities(
    # getTraceLength(

    # match_all = np.apply_along_axis(lambda a :np.isclose(mts_coord[:,0],a[0],atol=tol),1,ms2_coord)  # boolean arrays indicating the (mis)matches in mst; shape = [ms2,mts]
    # where_matches = np.apply_along_axis(lambda a : np.where(a),1,match_all)  # does not work because each ms2 has different amount of matches; shape=[ms2,matching mts (49,60,...)]

