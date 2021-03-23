"""MotifSet class and MotifSetIterator class definiton. 

A MotifSet object contains a set of DNA motifs whose occurrences will be
serached on the given genome variation graph.

This class is useful to manage the case in which the user wants to 
find the occurrences of more than one single motif.
"""

from grafimo.motif import Motif

from typing import List


class MotifSet(object):
    """
    This class represents a set of DNA motifs given as PWMs.

    A MotifSet object is an iterabole container for Motif objects. It is 
    particularly useful when searching with a single set of genomic 
    coordinates more than just a single motif.

    In the case a MEME file conatining data of more than one motif, a
    MotifSet object will contain the data for each of the available
    DNA motf, storing them as different Motif objects.

    ...

    Attributes
    ----------
    _motifs : list
        list of motifs conatined in the current object
    _motifs_num : int
        number of Motif objects contained in the current object

    Methods
    -------
    addMotif(mtfs : List[Motif])
        add the motifs in the list to the MotifSet
    getMotifsList()
        returns a list object containing all the motifs contained in
        the MotifSet
    length()
        returns the length of the MotifSet or the number of motifs 
        currently contained in the Motifset 
    """

    #-------------------------------------------------------------------
    # MotifSet attributes
    #
    _motifs: List[Motif] = list()
    _motifs_num: int = 0

    #-------------------------------------------------------------------
    # MotifSet methods
    #

    # these errors should never appear --> no need for error formatting
    # can assume that debug mode == True
    def __init__(self):
        pass


    def __iter__(self):
        return MotifSetIterator(self)

    
    def __len__(self):
        return len(self._motifs)


    def addMotif(self, mtfs: List[Motif]) -> None:
        """Add a list object containing Motif instanced to the Motifset
        object.

        ...

        Parameters
        ----------
        mtfs : list
            list containing of Motif objects
        """

        errmsg: str
        if not isinstance(mtfs, list):
            errmsg = "\n\nERROR: Expected list, got {}.\n"
            raise TypeError(errmsg.format(type(mtfs).__name__))
        if not all(isinstance(m, Motif) for m in mtfs):
            errmsg = "\n\nERROR: Expected Motif, got {}.\n"
            raise TypeError(errmsg.format(type(m).__name__))
        self._motifs += mtfs  # add motifs to self
        oldsize = self._motifs_num   
        self._motifs_num += len(mtfs)  # update the motif set size
        assert self._motifs_num > 0
        assert self._motifs_num == (oldsize + len(mtfs))
        

    def getMotifsList(self) -> List[Motif]:
        """Returns a list containing the Motif instances of MotifSet.

        ...

        Returns
        -------
        list
            list of Motif objects contained in MotifSet
        """

        errmsg: str
        if not self._motifs:
            errmsg = "\n\nERROR: The MotifSet object is empty.\n"
            raise ValueError(errmsg)
        if self._motifs_num == 0:
            errmsg = "\n\nERROR: The MotifSet object is empty.\n"
            raise ValueError(errmsg)
        return self._motifs
   

    def _size(self) -> int:
        """Compute the MotifSet object size.

        ...

        Returns
        -------
        int
            number of motifs in MotifSet
        """

        assert self._motifs_num >= 0
        return self._motifs_num

    @property
    def size(self):
        return self._size()

# end of MotifSet


class MotifSetIterator:
    """
    Iterator definition for MotifSet class. 

    ...

    Attributes
    ----------
    _motifSet : MotifSet
        reference to the itereted MotifSet
    _index : int
        element index
    """

    #-------------------------------------------------------------------
    # MotifSetIterator attributes
    #
    _motifSet: MotifSet
    _index: int

    #-------------------------------------------------------------------
    # MotifSetIterator methods
    #

    # these errors should never appear --> no need for error formatting
    # can assume that debug mode == True
    def __init__(self, motifSet: MotifSet):
        errmsg: str
        if not isinstance(motifSet, MotifSet):
            errmsg = "\n\nERROR: Expected MotifSet, got {}"
            raise ValueError(errmsg.format(type(motifSet).__name__))
        self._motifSet = motifSet
        self._index = 0


    def __next__(self) -> Motif:
        if self._index < self._motifSet.size:
            motifs_lst: List[Motif] = self._motifSet.getMotifsList()
            result: Motif = motifs_lst[self._index]
            self._index += 1
            return result
        raise StopIteration
   
# end of MotifSetIterator

