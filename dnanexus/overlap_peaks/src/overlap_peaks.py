#!/usr/bin/env python
# overlap_peaks 0.0.1

import re
import dxpy
import common
import logging

logger = logging.getLogger(__name__)
logger.addHandler(dxpy.DXLogHandler())
logger.propagate = False
logger.setLevel(logging.INFO)


@dxpy.entry_point('main')
def main(rep1_peaks, rep2_peaks, pooled_peaks,
         chrom_sizes, as_file, peak_type,
         pooledpr1_peaks=None, pooledpr2_peaks=None,
         prefix=None,
         rep1_signal=None, rep2_signal=None, pooled_signal=None):

    replicate_inputs = [pooledpr1_peaks, pooledpr2_peaks]
    simplicate_analysis = \
        all([replicate_input is None for replicate_input in replicate_inputs])

    rep1_peaks      = dxpy.DXFile(rep1_peaks)
    rep2_peaks      = dxpy.DXFile(rep2_peaks)
    pooled_peaks    = dxpy.DXFile(pooled_peaks)
    if not simplicate_analysis:
        pooledpr1_peaks = dxpy.DXFile(pooledpr1_peaks)
        pooledpr2_peaks = dxpy.DXFile(pooledpr2_peaks)
    chrom_sizes     = dxpy.DXFile(chrom_sizes)
    as_file         = dxpy.DXFile(as_file)

    # Input filenames - necessary to define each explicitly because input files
    # could have the same name, in which case subsequent
    # file would overwrite previous file
    rep1_peaks_fn      = 'rep1-%s' % (rep1_peaks.name)
    rep2_peaks_fn      = 'rep2-%s' % (rep2_peaks.name)
    pooled_peaks_fn    = 'pooled-%s' % (pooled_peaks.name)
    if not simplicate_analysis:
        pooledpr1_peaks_fn = 'pooledpr1-%s' % (pooledpr1_peaks.name)
        pooledpr2_peaks_fn = 'pooledpr2-%s' % (pooledpr2_peaks.name)
    chrom_sizes_fn     = 'chrom.sizes'
    as_file_fn         = '%s.as' % (peak_type)

    # Output filenames
    if prefix:
        basename = prefix
    else:
        # strip off the peak and compression extensions
        m = re.match(
            '(.*)(\.%s)+(\.((gz)|(Z)|(bz)|(bz2)))' % (peak_type),
            pooled_peaks.name)
        if m:
            basename = m.group(1)
        else:
            basename = pooled_peaks.name

    overlapping_peaks_fn    = '%s.replicated.%s' % (basename, peak_type)
    overlapping_peaks_bb_fn = overlapping_peaks_fn + '.bb'
    rejected_peaks_fn       = '%s.rejected.%s' % (basename, peak_type)
    rejected_peaks_bb_fn    = rejected_peaks_fn + '.bb'

    # Intermediate filenames
    overlap_tr_fn = 'replicated_tr.%s' % (peak_type)
    overlap_pr_fn = 'replicated_pr.%s' % (peak_type)

    # Download file inputs to the local file system with local filenames

    dxpy.download_dxfile(rep1_peaks.get_id(), rep1_peaks_fn)
    dxpy.download_dxfile(rep2_peaks.get_id(), rep2_peaks_fn)
    dxpy.download_dxfile(pooled_peaks.get_id(), pooled_peaks_fn)
    if simplicate_analysis:
        dxpy.download_dxfile(pooledpr1_peaks.get_id(), pooledpr1_peaks_fn)
        dxpy.download_dxfile(pooledpr2_peaks.get_id(), pooledpr2_peaks_fn)
    dxpy.download_dxfile(chrom_sizes.get_id(), chrom_sizes_fn)
    dxpy.download_dxfile(as_file.get_id(), as_file_fn)

    # the only difference between the peak_types is how the extra columns are
    # handled
    if peak_type == "narrowPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-10'
        bed_type = 'bed6+4'
    elif peak_type == "gappedPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$18-$17; if (($31/s1 >= 0.5) || ($31/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-15'
        bed_type = 'bed12+3'
    elif peak_type == "broadPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-9'
        bed_type = 'bed6+3'
    else:
        assert peak_type in ['narrowPeak', 'gappedPeak', 'broadPeak'], "%s is unrecognized.  peak_type should be narrowPeak, gappedPeak or broadPeak." % (peak_type)

    # Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as
    # the fractional overlap wrt any one of the overlapping peak pairs  > 0.5
    out, err = common.run_pipe([
        'intersectBed -wo -a %s -b %s' % (pooled_peaks_fn, rep1_peaks_fn),
        awk_command,
        cut_command,
        'sort -u',
        'intersectBed -wo -a stdin -b %s' % (rep2_peaks_fn),
        awk_command,
        cut_command,
        'sort -u'
        ], overlap_tr_fn)
    print(
        "%d peaks overlap with both true replicates"
        % (common.count_lines(overlap_tr_fn)))

    if not simplicate_analysis:
        # Find pooled peaks that overlap PseudoRep1 and PseudoRep2 where
        # overlap is defined as the fractional overlap wrt any one of the
        # overlapping peak pairs  > 0.5
        out, err = common.run_pipe([
            'intersectBed -wo -a %s -b %s' % (pooled_peaks_fn, pooledpr1_peaks_fn),
            awk_command,
            cut_command,
            'sort -u',
            'intersectBed -wo -a stdin -b %s' % (pooledpr2_peaks_fn),
            awk_command,
            cut_command,
            'sort -u'
            ], overlap_pr_fn)
        print(
            "%d peaks overlap with both pooled pseudoreplicates"
            % (common.count_lines(overlap_pr_fn)))

        # Combine peak lists
        out, err = common.run_pipe([
            'cat %s %s' % (overlap_tr_fn, overlap_pr_fn),
            'sort -u'
            ], overlapping_peaks_fn)
        print(
            "%d peaks overlap with true replicates or with pooled pseudoreplicates"
            % (common.count_lines(overlapping_peaks_fn)))
    else:  # this is a simplicate analysis
        # overlapping peaks are just based on pseudoreps of the one pool
        out, err = common.run_pipe([
            'cat %s' % (overlap_tr_fn),
            'sort -u'
            ], overlapping_peaks_fn)
        print(
            "%d peaks overlap"
            % (common.count_lines(overlapping_peaks_fn)))

    # rejected peaks
    out, err = common.run_pipe([
        'intersectBed -wa -v -a %s -b %s' % (pooled_peaks_fn, overlapping_peaks_fn)
        ], rejected_peaks_fn)
    print("%d peaks were rejected" % (common.count_lines(rejected_peaks_fn)))

    npeaks_in       = common.count_lines(common.uncompress(pooled_peaks_fn))
    npeaks_out      = common.count_lines(overlapping_peaks_fn)
    npeaks_rejected = common.count_lines(rejected_peaks_fn)

    # make bigBed files for visualization
    overlapping_peaks_bb_fn = common.bed2bb(
        overlapping_peaks_fn, chrom_sizes_fn, as_file_fn, bed_type=bed_type)
    rejected_peaks_bb_fn = common.bed2bb(
        rejected_peaks_fn, chrom_sizes_fn, as_file_fn, bed_type=bed_type)

    # Upload file outputs from the local file system.

    overlapping_peaks       = dxpy.upload_local_file(common.compress(overlapping_peaks_fn))
    overlapping_peaks_bb    = dxpy.upload_local_file(overlapping_peaks_bb_fn)
    rejected_peaks          = dxpy.upload_local_file(common.compress(rejected_peaks_fn))
    rejected_peaks_bb       = dxpy.upload_local_file(rejected_peaks_bb_fn)

    output = {
        "overlapping_peaks"     : dxpy.dxlink(overlapping_peaks),
        "overlapping_peaks_bb"  : dxpy.dxlink(overlapping_peaks_bb),
        "rejected_peaks"        : dxpy.dxlink(rejected_peaks),
        "rejected_peaks_bb"     : dxpy.dxlink(rejected_peaks_bb),
        "npeaks_in"             : npeaks_in,
        "npeaks_out"            : npeaks_out,
        'npeaks_rejected'       : npeaks_rejected
    }

    # These are just passed through for convenience so that signals and tracks
    # are available in one place.  Both input and output are optional.
    if rep1_signal:
        output.update({"rep1_signal": rep1_signal})
    if rep2_signal:
        output.update({"rep2_signal": rep2_signal})
    if pooled_signal:
        output.update({"pooled_signal": pooled_signal})

    return output


dxpy.run()
