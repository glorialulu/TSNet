=====
Get started
=====

To use wmoc in a project::

    import wmoc


Simple example
---------------

Simulate valve closure in a simple network.

.. code-block:: python

    # open an example network
    inp_file = 'LoopedNet_pump_valve.inp'

    # name of the valve to be close
    valve_to_close = ['11']
    # valve operation rule  
    tc = 2 # valve closure peroid
    ts = 0 # valve closure start time 
    se = 0 # end open percentage 
    m = 1 # closure constant
    valve_op = [tc,ts,se,m]

    wn, H1, V1, tt = wmoc.simulation.MOCSimulator(inp_file, dt, tf, 
                            valve_to_close, valve_op) 

    # specify the pipe name to visulize results
    pipe = '0'
    wmoc.postprocessing.plot_head_history(pipe,H1,wn,tt)
