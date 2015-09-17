
import constants

class lj_pot:
    def __init__(self):
        self.sigma = 2.315
        self.epsilon = 0.167
        self.mass = 63.55 * constants.amuToInternalMass
        self.lat = 3.615
        self.lattice_type = 'FCC'
        self.cutoff = 2.5*self.sigma
        self.name = "Cu"
        self.atomic_no = 29

    def output(self):
        pass

    def computeForce(self, atoms):
        POT_SHIFT = 1.0
        #POT_SHIFT = 0.0
        rCut2 = self.cutoff * self.cutoff

        for i in range(atoms.nAtoms):
            atoms.f[i,:] = 0.0

        ePot = 0.0

        s6 = self.sigma * self.sigma * self.sigma * self.sigma * self.sigma * self.sigma

        rCut6 = s6/(rCut2*rCut2*rCut2)
        eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0)

        # loop over atoms

        iPot = 0
        for i in range(atoms.nAtoms):
            for j in range(atoms.nAtoms):
                if i < j:
                    r2 = 0.0
                    dx = atoms.r[i,0] - atoms.r[j,0]
                    if dx < -0.5*atoms.bounds[0]:
                        dx += atoms.bounds[0]
                    if dx > 0.5*atoms.bounds[0]:
                        dx -= atoms.bounds[0]
                    r2 += dx*dx

                    dy = atoms.r[i,1] - atoms.r[j,1]
                    if dy < -0.5*atoms.bounds[1]:
                        dy += atoms.bounds[1]
                    if dy > 0.5*atoms.bounds[1]:
                        dy -= atoms.bounds[1]
                    r2 += dy*dy

                    dz = atoms.r[i,2] - atoms.r[j,2]
                    if dz < -0.5*atoms.bounds[2]:
                        dz += atoms.bounds[2]
                    if dz > 0.5*atoms.bounds[2]:
                        dz -= atoms.bounds[2]
                    r2 += dz*dz

                    if r2 > rCut2:
                        continue

                    ir2 = 1.0/r2

                    r6 = s6 * (ir2*ir2*ir2)
                    eLocal = r6*(r6-1.0) - eShift

                    ePot += eLocal
                    iPot += 1


                    fr = -4.0*self.epsilon * r6 * ir2*(12.0*r6 - 6.0)

                    atoms.f[i,0] -= dx*fr
                    atoms.f[j,0] += dx*fr
                    atoms.f[i,1] -= dy*fr
                    atoms.f[j,1] += dy*fr
                    atoms.f[i,2] -= dz*fr
                    atoms.f[j,2] += dz*fr

        ePot = ePot*4.0*self.epsilon
        return ePot

