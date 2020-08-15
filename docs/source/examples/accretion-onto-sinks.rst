--------------------
Accretion onto sinks
--------------------

Plot mass accretion and accretion rate onto sink particles.

.. figure:: ../_static/accretion.png

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    import plonk

    # Load simulation
    sim = plonk.load_sim(prefix='disc')
    sink_labels = ('Star', 'Planet')

    # Initialize figure
    fig, ax = plt.subplots(ncols=1, nrows=2, sharex=True, figsize=(12, 10))

    # Loop over sinks and plot
    for idx, sink in enumerate(sim.time_series['sinks']):

        # Time in years
        sink['time [year]'] = (sink['time [s]'].to_numpy() * plonk.units['s']).to('year').m

        # Mass accreted in Earth mass
        sink['mass_accreted [earth_mass]'] = (
            (sink['mass_accreted [kg]'].to_numpy() * plonk.units['kg'])
            .to('earth_mass')
            .magnitude
        )

        # Calculate accretion rate
        sink['accretion_rate [kg / s]'] = np.gradient(
            sink['mass_accreted [kg]'], sink['time [s]']
        )

        # Take rolling average
        accretion_rate = (
            sink['accretion_rate [kg / s]'].rolling(window=100).mean().to_numpy()
        )

        # Convert to Earth mass / year
        sink['accretion_rate [earth_mass / year]'] = (
            (accretion_rate * plonk.units['kg/s']).to('earth_mass / year').magnitude
        )

        # Plot
        sink.plot(
            'time [year]',
            'mass_accreted [earth_mass]',
            ax=ax[0],
            label=f'{sink_labels[idx]}',
        )
        sink.plot(
            'time [year]',
            'accretion_rate [earth_mass / year]',
            ax=ax[1],
            label=f'{sink_labels[idx]}',
        )

    # Set plot labels
    ax[0].set_ylabel('Mass accreted [$M_{\oplus}$]')
    ax[1].set_xlabel('Time [yr]')
    ax[1].set_ylabel('Accretion rate [$M_{\oplus}$/yr]')

    plt.show()
