"""
QA figure of the photometric calibration of a night.
"""

import numpy as np

from matplotlib.figure import Figure

from pipnick.qa.display import save_figure


def plot_residuals(meas, content, ofile):
    """
    Write the QA figure of the calibration of a night: the residuals of each
    calibrated band against airmass, color, catalog magnitude, and time.

    Used measurements are filled circles, rejected measurements are crosses,
    and unusable measurements are not shown.  Bands whose extinction was not
    measured are labeled.  Nothing is written if no band is calibrated.

    Parameters
    ----------
    meas : :class:`astropy.table.Table`
        The standard-star measurements; see
        :func:`~pipnick.pipelines.calibrate.calibrate_night`.
    content : dict
        The calibration; see :func:`~pipnick.pipelines.calibrate.write_calibration`.
    ofile : str, :class:`pathlib.Path`
        The output file.
    """
    bands = [b for b, e in content['band'].items() if e.get('status') == 'ok']
    if len(bands) == 0:
        return
    fig = Figure(figsize=(14, 3 * len(bands)))
    axes = fig.subplots(len(bands), 4, squeeze=False)
    mjd0 = np.min(meas['mjd']) if len(meas) > 0 else 0.
    for row, band in zip(axes, bands):
        entry = content['band'][band]
        sel = meas['band'] == band
        resid = np.ma.filled(meas['resid'], 0.)
        err = np.ma.filled(meas['resid_err'], 0.)
        xs = {
            'airmass': np.asarray(meas['airmass']),
            entry['color']: np.ma.filled(meas['color'], 0.),
            f'{band} (catalog)': np.ma.filled(meas['mag_cat'], 0.),
            'hours': (np.asarray(meas['mjd']) - mjd0) * 24.
        }
        used = sel & np.asarray(meas['used'])
        rejected = sel & np.asarray(meas['rejected'])
        for ax, (label, x) in zip(row, xs.items()):
            ax.errorbar(x[used], resid[used], yerr=err[used], fmt='o', color='C0', ms=4)
            ax.scatter(x[rejected], resid[rejected], marker='x', color='C3', s=30)
            ax.axhline(0., color='k', lw=0.5)
            ax.set_xlabel(label)
            ax.set_ylabel(f'{band} residual (mag)')
        title = f'{band}: rms {entry["rms"]:.3f} mag, {entry["n_stars"]} stars'
        if 'rms_loo' in entry:
            title += f', leave-one-out rms {entry["rms_loo"]:.3f}'
        row[0].set_title(title, loc='left', fontsize=9)
        if 'fallback' in entry['source']:
            for ax in row:
                ax.text(
                    0.5, 0.92, 'EXTINCTION NOT MEASURED: nominal value used',
                    transform=ax.transAxes, ha='center', va='top', color='C3', fontsize=8,
                    fontweight='bold'
                )
        if 'drift' in entry:
            row[3].set_title(
                f'drift {entry["drift"][0]:+.4f} +/- {entry["drift"][1]:.4f} mag/h',
                loc='right', fontsize=9
            )
    fig.suptitle(
        f'{content["calibration"]["setup"]} {content["calibration"]["night"]}: calibrated minus '
        'catalog magnitudes'
    )
    fig.tight_layout()
    save_figure(fig, ofile)
