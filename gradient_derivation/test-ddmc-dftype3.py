from ase.calculators.dftb import Dftb
from ase.calculators.ddmc import dDMC
from ase import Atoms
from ase.optimize import QuasiNewton
from ase.data.molecules import molecule
from ase.io import write
from ase.io import read


def numeric_force(atoms, a, i, d=0.001):
    """Evaluate force along i'th axis on a'th atom using finite difference.

    This will trigger two calls to get_potential_energy(), with atom a moved
    plus/minus d in the i'th axial direction, respectively.
    """
    p0 = atoms.positions[a, i]
    atoms.positions[a, i] += d
    eplus = atoms.get_potential_energy()
    print 
    atoms.positions[a, i] -= 2 * d
    eminus = atoms.get_potential_energy()
    atoms.positions[a, i] = p0
    return (eminus - eplus) / (2 * d)

def numeric_forces(atoms, indices=None, axes=(0, 1, 2), d=0.001,
                   parallel=None, name=None):
    print "BANANA1"
    """Evaluate finite-difference forces on several atoms.

    Returns an array of forces for each specified atomic index and
    each specified axis, calculated using finite difference on each
    atom and direction separately.  Array has same shape as if
    returned from atoms.get_forces(); uncalculated elements are zero.

    Calculates all forces by default."""
    
    import numpy as np
    from ase.parallel import world, rank, distribute_cpus
    from ase.utils import opencew


    if indices is None:
        indices = range(len(atoms))
    F_ai = np.zeros_like(atoms.positions)
    n = len(indices) * len(axes)
    if parallel is None:
        atom_tasks = [atoms] * n
        master = True
        calc_comm = world
    else:
        calc_comm, tasks_comm, tasks_rank = distribute_cpus(parallel, world)
        master = calc_comm.rank == 0
        calculator = atoms.get_calculator()
        calculator.set(communicator=calc_comm)
        atom_tasks = [None] * n
        for i in range(n):
            if ((i - tasks_rank) % tasks_comm.size) == 0:
                atom_tasks[i] = atoms
    for ia, a in enumerate(indices):
        for ii, i in enumerate(axes):
            atoms = atom_tasks[ia * len(axes) + ii]
            if atoms is not None:
                done = 0
                if name:
                    fname = '%s.%d%s.pckl' % (name, a, 'xyz'[i])
                    fd = opencew(fname, calc_comm)
                    if fd is None:
                        if master:
                            try:
                                F_ai[a, i] = pickle.load(open(fname))
                                print '# atom', a, 'xyz'[i], 'done'
                                done = 1
                            except EOFError:
                                pass
                        done = calc_comm.sum(done)
                if not done:
                    print '# rank', rank, 'calculating atom', a, 'xyz'[i]
                    force = numeric_force(atoms, a, i, d)
                    if master:
                        F_ai[a, i] = force
                        if name:
                            fd = open('%s.%d%s.pckl' % (name, a, 'xyz'[i]),
                                      'w')
                            pickle.dump(force, fd)
                            fd.close()
    if parallel is not None:
        world.sum(F_ai)
    return F_ai


#test = Atoms('2N',positions=[(0.,0.,0.),(0.,0.,4.10)])
#test = molecule('C2H2')
test = read('S668-24BenzeneBenzenepipi-090.xyz')
#test = Atoms('2Ar', position=[(0.,0.,0.),(0.,0.,3.)])

    
test.set_calculator(dDMC(label='S668-24BenzeneBenzenepipi-090',atoms=test,dftype='3',tagtype='column'))

print test.get_potential_energy()
numer =  numeric_forces(test, indices=None, axes=(0, 1, 2), d=0.001,
                   parallel=None, name=None)

anal = test.get_forces()

print 'Numerical Gradient: '
print numer
print
print 'Analitical Gradient: '
print anal
print
print 'Numerical/Analytical:'
print numer/anal


