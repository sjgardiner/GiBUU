#!/usr/bin/env python

"""
Parser for the LHEF format

For more detail about the LHE Format:

[1] "Generic User Process Interface for Event Generators"
    http://arxiv.org/abs/hep-ph/0109068

    Original Les Houches Accord for user-defined processes, using Fortran
    commonblocks.

[2] "A Standard Format for Les Houches Events"
    http://arxiv.org/abs/hep-ph/0609017

    The first proposal for the .lhe format


"""

from textwrap import dedent
from pdg_codes import PDG_CODES

__author__ = "Roberto Vidal"
__copyright__ = "(c) 2014 by José Roberto Vidal"
__version__ = "0.0.2"
__license__ = "GPL"


def iterate(lines):
    """
    Reads and parses a LHEF file and returns an iterator for the events
    contained therein.

    lines can be an iterable containing the lines or a filename
    The iterator gets exhausted after one full iteration!
    """
    is_file = isinstance(lines, str)
    if is_file:
        lines = open(lines)

    buff = []
    text_buff = []
    buffering = False
    try:
        for (lineno, line) in enumerate(lines):
            line = line.strip()

            if line.find(_CLOSE_EVENT_TAG) > -1:
                yield Event(buff)

                buff = []
                text_buff = []
                buffering = False

            if buffering:
                line = line.strip()
                text_buff.append(line)
                buff.append(line)

            if line.find(_OPEN_EVENT_TAG) > -1:
                buffering = True

    except _EventParseError as e:
        err_line = lineno - len(buff) + e.index + 1
        ex = ParseError(e.message, lineno=err_line, text=text_buff[e.index])
        raise ex from None
    finally:
        if is_file:
            lines.close()


def read(lines):
    """
    Reads and parses a LHEF file and returns a LhefFile object with all the
    contents loaded in memory

    lines can be an iterable containing the lines or a filename
    """
    leh_file = LhefFile()
    try:
        leh_file.events = list(iterate(lines))
    except ParseError as e:
        if (isinstance(lines, str)):
            e.filename = lines
        raise
    else:
        return leh_file


class LhefFile:
    """
    Represents a parsed .lhe file
    Exposed fields:
        events[]
    """
    def __init__(self, raw_events=None):
        if raw_events is not None:
            self.events = list(map(Event, raw_events))

    def __repr__(self):
        return "<LhefFile, events: %d>" % len(self.events)


class Event:
    """
    Represents an event
    Exposed fields:
        particles[], nup, process_id, weight, scale, alpha_qcd, alpha_qed
        opt_info[]

    opt_info stores the textual optional lines, since their format is
    completely unspecified [2]
    """
    def __init__(self, event_lines):
        """
        Constructs an event object from the textual version event_lines
        """
        metadata = list(map(float, (event_lines[0].split())))
        nup = int(metadata[0])

        if nup > len(event_lines) - 1:
            raise _EventParseError("Inconsistent number of particles", 0)

        if len(metadata) != Event._FIELDS:
            raise _EventParseError(
                "Encountered %d event fields instead of %d" %
                (len(metadata), Event._FIELDS), 0
            )

        particles = map(Event._parse_event_line, event_lines[1:nup + 1])
        self.particles = [None] * nup

        try:
            for (i, raw_particle) in enumerate(particles):
                self.particles[i] = Particle(raw_particle)
        except _EventParseError as e:
            e.index = i
            raise

        if nup < len(event_lines) - 1:
            self.opt_info = event_lines[nup + 1:]

        # NUP
        self.nup = nup
        # IDPRUP
        self.process_id = int(metadata[1])
        # XWGTUP
        self.weight = metadata[2]
        # SCALUP scale of the event in GeV
        self.scale = metadata[3]
        # AQCDUP
        self.alpha_qcd = metadata[4]
        # AQEDUP
        self.alpha_qed = metadata[5]

    def __repr__(self):
        return "<Event, particles: %d>" % len(self.particles)

    _FIELDS = 6

    def _parse_event_line(s):
        return list(map(float, s.split()))


class Particle:
    """
    Represents a particle within an event
    Exposed fields:
        id, status, first_mother, last_mother, color(), p(), mass,
        lifetime, spin
    """
    def __init__(self, raw_particle):

        if len(raw_particle) != Particle._FIELDS:
            raise _EventParseError(
                "Encountered %d fields for a particle instead of %d" %
                (len(raw_particle), Particle._FIELDS)
            )

        # IDUP
        self.id = int(raw_particle[0])
        self._name = Particle._name(self.id)
        # ISTUP
        self.status = int(raw_particle[1])

        # MOTHUP(1)
        self.first_mother = int(raw_particle[2])

        # MOTHUP(2)
        # [1] index of last mother (for decays, particles
        # have one mother only, so this is zero)
        self.last_mother = int(raw_particle[3])

        # ICOLUP 1, 2
        self.color = (int(raw_particle[4]), int(raw_particle[5]))

        # PUP(4, 1-3)
        # 4-momentum (E, px, py, pz)
        self.p = (
            raw_particle[9], raw_particle[6], raw_particle[7], raw_particle[8]
        )

        # PUP(5)
        self.mass = raw_particle[10]

        # VTIMUP
        # [1] invariant lifetime c\tau (distance from production
        # to decay) in mm
        self.lifetime = raw_particle[11]

        # SPINUP
        # [1] cosine of the angle between the spin-vector of particle I and the
        # 3-momentum of the decaying particle, specified in the lab frame
        self.spin = raw_particle[12]

    def __repr__(self):
        return "<Particle %s, mass: %1.e, momentum: %s>" % (
            self._name, self.mass, "(%1.e,%1.e,%1.e,%1.e)" % self.p
        )

    def __str__(self):
        return dedent("""
            %(name)s
                * Energy: %(energy)f GeV
                * Mass: %(mass)f GeV
                * Trimomentum: %(momentum)s GeV
                * Color: %(color)s
                * Status: %(status)s
                * Mothers: %(first_mother)d, %(last_mother)d
                * Spin: %(spin)s
                * Lifetime: %(lifetime)f mm
                * PDG number: %(id)d
        """).strip() % {
            "id": self.id,
            "mass": self.mass,
            "momentum": "(%f, %f, %f)" % self.p[1:],
            "energy": self.p[0],
            "color": "(%d, %d)" % self.color,
            "status": _STATUS_NAMES[self.status],
            "first_mother": self.first_mother,
            "last_mother": self.last_mother,
            "spin": "unknown" if self.spin == 9 else "%f" % self.spin,
            "lifetime": self.lifetime,
            "name": self._name
        }

    _FIELDS = 13

    # Helper to get the particle name
    def _name(id):
        if id in PDG_CODES:
            return PDG_CODES[id]
        else:
            return "<unknown %d>" % id


class ParseError(Exception):
    """
    Exception to be raised if parsing fails
    """
    def __init__(self, message, *, lineno=None, text=None, filename=None):
        self.message = message
        self.lineno = lineno
        self.text = text
        self.filename = filename

    def __str__(self):
        return dedent("""
            %s in line %d%s:

            > %s
        """).rstrip() % (
            self.message, self.lineno,
            " of file %s" % self.filename if self.filename is not None else "",
            self.text,
        )


class _EventParseError(ParseError):
    def __init__(self, message, index=None):
        self.message = message
        self.index = index

    def __str__(self):
        pass


# These are the possible values of particle.status
# Accessible as Particle.INCOMING_PARTICLE, etc.
INCOMING_PARTICLE = -1
OUTGOING_PARTICLE = +1
INTERMEDIATE_SPACELIKE = -2
INTERMEDIATE_RESONANCE = +2
INTERMEDIATE_RESONANCE_DOC_ONLY = +3
INCOMING_BEAM = -9


_STATUS_NAMES = {
    -1: "INCOMING_PARTICLE",
    +1: "OUTGOING_PARTICLE",
    -2: "INTERMEDIATE_SPACELIKE",
    +2: "INTERMEDIATE_RESONANCE",
    +3: "INTERMEDIATE_RESONANCE_DOC_ONLY",
    -9: "INCOMING_BEAM",
}

_OPEN_EVENT_TAG = "<event"
_CLOSE_EVENT_TAG = "</event"


__all__ = [
    "read", "Event", "Particle", "LhefFile", "iterate", "ParseError",
    "INCOMING_PARTICLE",
    "OUTGOING_PARTICLE",
    "INTERMEDIATE_SPACELIKE",
    "INTERMEDIATE_RESONANCE",
    "INTERMEDIATE_RESONANCE_DOC_ONLY",
    "INCOMING_BEAM",
]
