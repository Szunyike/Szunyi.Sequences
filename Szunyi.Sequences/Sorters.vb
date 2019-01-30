Imports Bio
Imports Bio.IO.GenBank


Public Class Sorters
    Public Class ByID
        Implements IComparer(Of Bio.ISequence)

        Public Function Compare(x As ISequence, y As ISequence) As Integer Implements IComparer(Of ISequence).Compare
            Return x.ID.CompareTo(y.ID)
        End Function
    End Class
    Public Class BySeq
        Implements IComparer(Of Bio.ISequence)

        Public Function Compare(x As ISequence, y As ISequence) As Integer Implements IComparer(Of ISequence).Compare
            Dim min = System.Math.Min(x.Count, y.Count)
            For i1 = 0 To min - 1
                If x(i1) <> y(i1) Then Return x(i1).CompareTo(y(i1))
            Next
            Return x.Count.CompareTo(y.Count)
        End Function
    End Class
End Class
