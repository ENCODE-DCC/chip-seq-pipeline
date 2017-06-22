#!/usr/bin/env python2
# encode_spp 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# See https://wiki.dnanexus.com/Developer-Portal for documentation and
# tutorials on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import subprocess
import common
import dxpy
import logging

logger = logging.getLogger(__name__)
logger.addHandler(dxpy.DXLogHandler())
logger.propagate = False
logger.setLevel(logging.INFO)


def spp(experiment, control, xcor_scores, chrom_sizes, spp_version,
        bigbed=False, as_file=None, name="spp", prefix=None,
        fragment_length=None):
    spp_applet = \
        dxpy.find_one_data_object(
            classname='applet',
            name='spp',
            project=dxpy.PROJECT_CONTEXT_ID,
            zero_ok=False,
            more_ok=False,
            return_handler=True)
    spp_input = {"experiment": experiment,
                 "control": control,
                 "xcor_scores_input": xcor_scores,
                 "bigbed": bigbed,
                 "chrom_sizes": chrom_sizes,
                 "spp_version": spp_version}
    if fragment_length is not None:
        spp_input.update({"fragment_length": fragment_length})
    if bigbed and as_file:
        spp_input.update({"as_file": as_file})
    if prefix:
        spp_input.update({"prefix": prefix})
    return spp_applet.run(spp_input, name=name)


def xcor_only(tags, paired_end, spp_version, name='xcor_only'):
    xcor_only_applet = \
        dxpy.find_one_data_object(
            classname='applet',
            name='xcor_only',
            project=dxpy.PROJECT_CONTEXT_ID,
            zero_ok=False,
            more_ok=False,
            return_handler=True)
    return xcor_only_applet.run(
        {"input_tagAlign": tags,
         "paired_end": paired_end,
         "spp_version": spp_version},
        name=name)


@dxpy.entry_point('main')
def main(rep1_ta, ctl1_ta, rep1_xcor, rep1_paired_end,
         npeaks, nodups,  chrom_sizes, spp_version,
         rep2_ta=None, ctl2_ta=None, rep2_xcor=None, rep2_paired_end=None,
         as_file=None, idr_peaks=False, fragment_length=None):

    rep1_ta_file = dxpy.DXFile(rep1_ta)
    dxpy.download_dxfile(rep1_ta_file.get_id(), rep1_ta_file.name)
    rep1_ta_filename = rep1_ta_file.name
    ntags_rep1 = common.count_lines(rep1_ta_filename)

    simplicate_experiment = rep1_ta and not rep2_ta
    if simplicate_experiment:
        logger.info("No rep2 tags specified so processing as a simplicate experiment.")
    else:
        logger.info("Rep1 and rep2 tags specified so processing as a replicated experiment.")

    if not simplicate_experiment:
        assert rep1_paired_end == rep2_paired_end, 'Mixed PE/SE not supported'
        rep2_ta_file = dxpy.DXFile(rep2_ta)
        dxpy.download_dxfile(rep2_ta_file.get_id(), rep2_ta_file.name)
        rep2_ta_filename = rep2_ta_file.name
        ntags_rep2 = common.count_lines(rep2_ta_filename)
    paired_end = rep1_paired_end

    unary_control = (ctl1_ta == ctl2_ta) or not ctl2_ta
    ctl1_ta_file = dxpy.DXFile(ctl1_ta)
    dxpy.download_dxfile(ctl1_ta_file.get_id(), ctl1_ta_file.name)
    ctl1_ta_filename = ctl1_ta_file.name

    if not unary_control:
        ctl2_ta_file = dxpy.DXFile(ctl2_ta)
        dxpy.download_dxfile(ctl2_ta_file.get_id(), ctl2_ta_file.name)
        ctl2_ta_filename = ctl2_ta_file.name
    else:
        ctl2_ta_file = ctl1_ta_file
        ctl2_ta_filename = ctl1_ta_file.name

    ntags_ctl1 = common.count_lines(ctl1_ta_filename)
    ntags_ctl2 = common.count_lines(ctl2_ta_filename)
    rep1_control = ctl1_ta  # default.  May be changed later.
    rep1_ctl_msg = "control rep1"
    rep2_control = ctl2_ta  # default.  May be changed later.
    rep2_ctl_msg = "control rep2"

    rep_info = [(ntags_rep1, 'replicate 1', rep1_ta_filename)]
    if not simplicate_experiment:
        rep_info.append((ntags_rep2, 'replicate 2', rep2_ta_filename))
    rep_info.extend(
        [(ntags_ctl1, 'control 1', ctl1_ta_filename),
         (ntags_ctl2, 'control 2', ctl2_ta_filename)])
    for n, name, filename in rep_info:
        logger.info("Found %d tags in %s file %s" % (n, name, filename))

    subprocess.check_output('ls -l', shell=True, stderr=subprocess.STDOUT)

    if not simplicate_experiment:
        pool_applet = dxpy.find_one_data_object(
                classname='applet',
                name='pool',
                project=dxpy.PROJECT_CONTEXT_ID,
                zero_ok=False,
                more_ok=False,
                return_handler=True)
        pool_replicates_subjob = \
            pool_applet.run(
                {"inputs": [rep1_ta, rep2_ta],
                 "prefix": 'pooled_reps'},
                name='Pool replicates')
        pooled_replicates = pool_replicates_subjob.get_output_ref("pooled")
        pooled_replicates_xcor_subjob = \
            xcor_only(
                pooled_replicates,
                paired_end,
                spp_version,
                name='Pool cross-correlation')

    if unary_control:
        logger.info("Only one control supplied.")
        if not simplicate_experiment:
            logger.info("Using one control for both replicate 1 and 2 and for the pool.")
        control_for_pool = rep1_control
        pool_ctl_msg = "one control"
    else:
        pool_controls_subjob = pool_applet.run(
            {"inputs": [ctl1_ta, ctl2_ta],
             "prefix": "PL_ctls"},
            name='Pool controls')
        pooled_controls = pool_controls_subjob.get_output_ref("pooled")
        # always use the pooled controls for the pool
        control_for_pool = pooled_controls
        pool_ctl_msg = "pooled controls"

        # use the pooled controls for the reps depending on the ratio of rep to
        # control reads
        ratio_ctl_reads = float(ntags_ctl1)/float(ntags_ctl2)
        if ratio_ctl_reads < 1:
                ratio_ctl_reads = 1/ratio_ctl_reads
        ratio_cutoff = 1.2
        if ratio_ctl_reads > ratio_cutoff:
                logger.info(
                    "Number of reads in controls differ by > factor of %f. Using pooled controls."
                    % (ratio_cutoff))
                rep1_control = pooled_controls
                rep2_control = pooled_controls
        else:
                if ntags_ctl1 < ntags_rep1:
                        logger.info("Fewer reads in control replicate 1 than experiment replicate 1.  Using pooled controls for replicate 1.")
                        rep1_control = pooled_controls
                        rep1_ctl_msg = "pooled controls"
                elif not simplicate_experiment and ntags_ctl2 < ntags_rep2:
                        logger.info("Fewer reads in control replicate 2 than experiment replicate 2.  Using pooled controls for replicate 2.")
                        rep2_control = pooled_controls
                        rep2_ctl_msg = "pooled controls"
                else:
                    logger.info(
                        "Using distinct controls for replicate 1 and 2.")
                    rep1_control = ctl1_ta  # default.  May be changed later.
                    rep2_control = ctl2_ta  # default.  May be changed later.
                    rep1_ctl_msg = "control rep1"
                    rep2_ctl_msg = "control rep2"

    common_args = {
        'chrom_sizes': chrom_sizes,
        'spp_version': spp_version,
        'as_file':     as_file
        }
    if fragment_length is not None:
        common_args.update({'fragment_length': fragment_length})
    rep1_peaks_subjob = spp(
        rep1_ta,
        rep1_control,
        rep1_xcor,
        bigbed=True,
        name='Rep1 peaks vs %s' % (rep1_ctl_msg),
        prefix='R1', **common_args)

    if not simplicate_experiment:
        rep2_peaks_subjob = spp(
            rep2_ta,
            rep2_control,
            rep2_xcor,
            bigbed=True,
            name='Rep2 peaks vs %s' % (rep2_ctl_msg),
            prefix='R2', **common_args)

        pooled_peaks_subjob = spp(
            pooled_replicates,
            control_for_pool,
            pooled_replicates_xcor_subjob.get_output_ref("CC_scores_file"),
            bigbed=True,
            name='Pooled peaks vs %s' % (pool_ctl_msg),
            prefix='PL', **common_args)

    output = {
        'rep1_peaks':       rep1_peaks_subjob.get_output_ref("peaks"),
        'rep1_peaks_bb':    rep1_peaks_subjob.get_output_ref("peaks_bb"),
        'rep1_xcor_plot':   rep1_peaks_subjob.get_output_ref("xcor_plot"),
        'rep1_xcor_scores': rep1_peaks_subjob.get_output_ref("xcor_scores")
    }

    if not simplicate_experiment:
        output.update({
            'rep2_peaks':       rep2_peaks_subjob.get_output_ref("peaks"),
            'rep2_peaks_bb':    rep2_peaks_subjob.get_output_ref("peaks_bb"),
            'rep2_xcor_plot':   rep2_peaks_subjob.get_output_ref("xcor_plot"),
            'rep2_xcor_scores': rep2_peaks_subjob.get_output_ref("xcor_scores"),

            'pooled_peaks':       pooled_peaks_subjob.get_output_ref("peaks"),
            'pooled_peaks_bb':    pooled_peaks_subjob.get_output_ref("peaks_bb"),
            'pooled_xcor_plot':   pooled_peaks_subjob.get_output_ref("xcor_plot"),
            'pooled_xcor_scores': pooled_peaks_subjob.get_output_ref("xcor_scores")
        })

    if idr_peaks:  # also call peaks on pseudoreplicates for IDR
        pseudoreplicator_applet = \
            dxpy.find_one_data_object(
               classname='applet',
               name='pseudoreplicator',
               project=dxpy.PROJECT_CONTEXT_ID,
               zero_ok=False,
               more_ok=False,
               return_handler=True)

        rep1_pr_subjob = \
            pseudoreplicator_applet.run(
                {"input_tags": rep1_ta,
                 "prefix": 'R1PR'},
                name='Pseudoreplicate rep1 -> R1PR1,2')

        rep1pr1_peaks_subjob = spp(
            rep1_pr_subjob.get_output_ref("pseudoreplicate1"),
            rep1_control,
            rep1_xcor,
            bigbed=False,
            name='R1PR1 peaks vs %s' % (rep1_ctl_msg),
            prefix='R1PR1', **common_args)

        rep1pr2_peaks_subjob = spp(
            rep1_pr_subjob.get_output_ref("pseudoreplicate2"),
            rep1_control,
            rep1_xcor,
            bigbed=False,
            name='R1PR2 peaks vs %s' % (rep1_ctl_msg),
            prefix='R1PR2', **common_args)

        output.update({
            'rep1pr1_peaks':         rep1pr1_peaks_subjob.get_output_ref("peaks"),
            'rep1pr2_peaks':         rep1pr2_peaks_subjob.get_output_ref("peaks")
            })

        if not simplicate_experiment:
            rep2_pr_subjob = \
                pseudoreplicator_applet.run(
                    {"input_tags": rep2_ta,
                     "prefix": 'R2PR'},
                    name='Pseudoreplicate rep2 -> R2PR1,2')

            pool_pr1_subjob = pool_applet.run({
                "inputs": [
                    rep1_pr_subjob.get_output_ref("pseudoreplicate1"),
                    rep2_pr_subjob.get_output_ref("pseudoreplicate1")],
                "prefix": 'PPR1'},
                name='Pool R1PR1+R2PR1 -> PPR1')

            pool_pr2_subjob = pool_applet.run({
                "inputs": [
                    rep1_pr_subjob.get_output_ref("pseudoreplicate2"),
                    rep2_pr_subjob.get_output_ref("pseudoreplicate2")],
                "prefix": 'PPR2'},
                name='Pool R1PR2+R2PR2 -> PPR2')

            rep2pr1_peaks_subjob = spp(
                rep2_pr_subjob.get_output_ref("pseudoreplicate1"),
                rep2_control,
                rep2_xcor,
                bigbed=False,
                name='R2PR1 peaks vs %s' % (rep2_ctl_msg),
                prefix='R2PR1', **common_args)

            rep2pr2_peaks_subjob = spp(
                rep2_pr_subjob.get_output_ref("pseudoreplicate2"),
                rep2_control,
                rep2_xcor,
                bigbed=False,
                name='R2PR2 peaks vs %s' % (rep2_ctl_msg),
                prefix='R2PR2', **common_args)

            pooledpr1_peaks_subjob = spp(
                pool_pr1_subjob.get_output_ref("pooled"),
                control_for_pool,
                pooled_replicates_xcor_subjob.get_output_ref("CC_scores_file"),
                bigbed=False,
                name='PPR1 peaks vs %s' % (pool_ctl_msg),
                prefix='PPR1', **common_args)

            pooledpr2_peaks_subjob = spp(
                pool_pr2_subjob.get_output_ref("pooled"),
                control_for_pool,
                pooled_replicates_xcor_subjob.get_output_ref("CC_scores_file"),
                bigbed=False,
                name='PPR2 peaks vs %s' % (pool_ctl_msg),
                prefix='PPR2', **common_args)

            output.update({
                'rep2pr1_peaks':         rep2pr1_peaks_subjob.get_output_ref("peaks"),
                'rep2pr2_peaks':         rep2pr2_peaks_subjob.get_output_ref("peaks"),
                'pooledpr1_peaks':       pooledpr1_peaks_subjob.get_output_ref("peaks"),
                'pooledpr2_peaks':       pooledpr2_peaks_subjob.get_output_ref("peaks"),
            })

    return output

dxpy.run()
