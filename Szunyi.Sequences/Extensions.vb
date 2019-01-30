Imports System.Runtime.CompilerServices
Imports Bio
Imports Bio.IO.GenBank
Imports Szunyi.Common.Extensions
Imports Szunyi.Sequences

Public Module Extensions

    Public SequnceID_Sorter As New Szunyi.Sequences.Sorters.ByID
    Public SequnceSeq_Sorter As New Szunyi.Sequences.Sorters.BySeq
#Region "Misc"
    <Extension()>
    Public Function LongestRepetition(Seq As Bio.ISequence, Nucleotide As Byte) As Integer
        Dim Started As Boolean = False
        Dim cStart As Int16 = -1
        Dim MaxLength As Integer = -1
        For i1 = 0 To Seq.Count - 1
            If Seq(i1) = Nucleotide And Started = False Then
                cStart = i1
                Started = True
            ElseIf Seq(i1) = Nucleotide Then

            ElseIf Started = True Then ' 
                Dim l = i1 - cStart
                If l > MaxLength Then MaxLength = l
                Started = False
            End If
        Next
        If Started = True Then
            Dim l = Seq.Count - cStart
            If l > MaxLength Then MaxLength = l
        End If
        Return MaxLength
    End Function

    <Extension()>
    Public Function LongestRepetition(SeqS As IEnumerable(Of Bio.ISequence), Nucleotide As Byte) As Integer
        Dim MaxLength = -1
        For Each Seq In SeqS
            Dim l = Seq.LongestRepetition(Nucleotide)
            If l > MaxLength Then MaxLength = l
        Next
        Return MaxLength
    End Function

    <Extension()>
    Public Function ToPartitionTable(Seqs As List(Of Bio.ISequence), Type As String) As String
        Dim str As New System.Text.StringBuilder
        Dim cLength As Integer = 1
        For i1 = 0 To Seqs.Count - 1
            Dim Seq = Seqs(i1)
            str.Append(Type).Append(", gene").Append(i1).Append(" = ").Append(cLength).Append("-")
            cLength += Seq.Count
            str.Append(cLength).AppendLine()
        Next
        str.Length -= 2
        Return str.ToString
    End Function

#End Region

#Region "Concatenate"
    <Extension()>
    Public Function Concatenate(SeqA As Bio.Sequence, SeqB As Bio.Sequence) As Bio.Sequence
        Dim Alp = SeqA.Alphabet
        Dim b(SeqA.Count + SeqB.Count - 2) As Byte
        SeqA.CopyTo(b, 0, SeqA.Count)
        SeqB.CopyTo(b, SeqB.Count, SeqB.Count)
        Dim Seq As New Bio.Sequence(SeqA.Alphabet, b)
        Seq.ID = Seq.ID
        Seq.Metadata = SeqA.Metadata
        Return Seq

    End Function
    <Extension()>
    Public Function Concatenate(SeqA As Bio.ISequence, SeqB As Bio.ISequence) As Bio.Sequence

        Dim nSeq(SeqA.Count + SeqB.Count - 1) As Byte
        Dim A = SeqA.ToArray
        Dim B = SeqB.ToArray
        A.CopyTo(nSeq, 0)
        B.CopyTo(nSeq, SeqA.Count)
        Dim Seq As New Bio.Sequence(SeqA.Alphabet, nSeq)
        Seq.ID = Seq.ID
        Seq.Metadata = SeqA.Metadata
        Return Seq

    End Function
    <Extension()>
    Public Function Concatenate(Seqs As List(Of Bio.ISequence)) As Bio.Sequence

        Dim Counts = (From x In Seqs Select x.Count).Sum
        Dim nSeq(Counts - 1) As Byte
        Dim cIndex As Integer = 0
        For Each Seqi In Seqs
            Dim A = Seqi.ToArray
            A.CopyTo(nSeq, cIndex)
            cIndex += Seqi.Count
        Next

        Dim Seq As New Bio.Sequence(Bio.Alphabets.AmbiguousProtein, nSeq)

        Seq.ID = Seqs.First.ID

        Return Seq

    End Function
    <Extension()>
    Public Function Concatenate(Seqs As List(Of List(Of Bio.ISequence))) As List(Of Bio.ISequence)
        ' 1 Is Everything the same length
        ' 2 
        Dim dist_nofItemCounts = (From x In Seqs Select x.Count).Distinct
        Dim out As New List(Of Bio.ISequence)
        If dist_nofItemCounts.Count = 1 Then
            For i1 = 0 To Seqs.First.Count - 1
                Dim nSeqs As New List(Of Bio.ISequence)
                For i2 = 0 To Seqs.Count - 1
                    nSeqs.Add(Seqs(i2)(i1))
                Next
                out.Add(nSeqs.Concatenate)
            Next

        End If
        Return out
    End Function
#End Region

#Region "Translate"
    <Extension()>
    Public Function TranslateFull(Seqs As List(Of Bio.ISequence), Optional RemoveTerminalStop As Boolean = False) As List(Of Bio.ISequence)

        Dim out As New List(Of Bio.ISequence)
        For Each Seq In Seqs
            out.Add(TranslateFull(Seq, RemoveTerminalStop))
        Next
        Return out

    End Function
    <Extension()>
    Public Function TranslateFull(Seq As Bio.ISequence, Optional RemoveTerminalStop As Boolean = False) As Bio.Sequence
        Dim RNASeq = Bio.Algorithms.Translation.Transcription.Transcribe(Seq)
        Dim AASeq As Bio.ISequence = Bio.Algorithms.Translation.ProteinTranslation.Translate(RNASeq, 0)
        AASeq.ID = Seq.ID
        If RemoveTerminalStop = False Then
            Return AASeq
        Else
            Return AASeq.RemoveTerminalStop
        End If
    End Function

    <Extension()>
    Public Function RemoveTerminalStop(Seq As Bio.ISequence) As Bio.Sequence
        If Seq.Last <> Bio.Alphabets.AmbiguousProtein.Ter Then Return Seq
        Dim S As Bio.Sequence = Seq
        Dim NewSeq(S.Count - 2) As Byte
        S.CopyTo(NewSeq, 0, Seq.Count - 1)
        Dim a As New Bio.Sequence(Bio.Alphabets.AmbiguousProtein, NewSeq)
        a.ID = Seq.ID
        Return a
    End Function

#End Region

#Region "Accession, SOurce, TaxID etc ..."
    ''' <summary>
    ''' Return Primary Accession or String.Empty
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function Accesion(Seq As Bio.ISequence) As String
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Accession) = True OrElse md.Accession.Primary = True Then Return String.Empty
        Return md.Accession.Primary
    End Function
    ''' <summary>
    ''' Retrun Source.Organism.Genus + Source.Organism.Species
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function CommonName(Seq As Bio.ISequence) As String
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Accession) = True OrElse md.Accession.Primary = True Then Return String.Empty
        Return md.Accession.Primary
    End Function
    ''' <summary>
    ''' Retrun Source or Nothing
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function Source(Seq As Bio.ISequence) As Bio.IO.GenBank.SequenceSource
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Source) = True Then Return Nothing
        Return md.Source
    End Function
    ''' <summary>
    ''' Retrun TaxId or String.Empty
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function TaxId(Seq As Bio.ISequence) As String
        Dim Source = Seq.FeatureSource
        If IsNothing(Source) = True Then Return String.Empty
        Dim db = Source.Qualifiers("db_xref")
        Dim res = From x In db Where x.Contains("taxon")
        If res.Count > 0 Then
            Dim s1 = Split(res.First, "taxon:").Last
            Return s1.Trim(Chr(34))
        End If

        Return String.Empty
    End Function
    ''' <summary>
    ''' Return Strain or String.Empty
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function Strain(Seq As Bio.ISequence) As String
        Dim Source = Seq.FeatureSource
        If IsNothing(Source) = True Then Return String.Empty
        If Source.Qualifiers.ContainsKey("strain") Then
            Dim db = Source.Qualifiers("strain")
            Return db.First.aZ09_
        End If
        Return String.Empty
    End Function
    ''' <summary>
    ''' Return First Source FeatureItem or Nothing
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function FeatureSource(Seq As Bio.ISequence) As Bio.IO.GenBank.FeatureItem
        Dim md As Bio.IO.GenBank.GenBankMetadata = Seq.Metadata(Bio.Util.Helper.GenBankMetadataKey)
        If IsNothing(md) = True Then Return Nothing
        Dim source = From x In md.Features.All Where x.Key = "source"
        If source.Count > 0 Then
            Return source.First
        End If
        Return Nothing
    End Function
#End Region

#Region "GenBankMetaData"
    <Extension()>
    Public Sub GenBankMetaData_NACreate(ByRef Seq As Bio.ISequence)
        Dim x As New Bio.IO.GenBank.GenBankMetadata
        x.Locus = New Bio.IO.GenBank.GenBankLocusInfo()
        x.Locus.Date = Now
        x.Locus.MoleculeType = Seq.MoleculeType
        x.Locus.Name = Seq.ID
        x.Locus.SequenceLength = Seq.Count
        x.Locus.StrandTopology = Bio.IO.GenBank.SequenceStrandTopology.Linear
        x.Accession = New GenBankAccession()
        x.Accession.Primary = Seq.ID
        x.Source = New Bio.IO.GenBank.SequenceSource
        x.Source.CommonName = "Unknown."
        x.Source.Organism = New Bio.IO.GenBank.OrganismInfo
        x.Source.Organism.Species = "Unknown"
        x.Source.Organism.Genus = "Unknown."
        x.Features = New SequenceFeatures
        Seq.Metadata.Add(Bio.Util.Helper.GenBankMetadataKey, x)
    End Sub
    <Extension>
    Private Iterator Function GenBankMetaDatas(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of Bio.IO.GenBank.GenBankMetadata)
        Dim Out As New List(Of Bio.IO.GenBank.GenBankMetadata)
        For Each Seq In Seqs
            If Seq.Metadata.ContainsKey(Bio.Util.Helper.GenBankMetadataKey) Then
                Yield (Seq.Metadata(Bio.Util.Helper.GenBankMetadataKey))
            End If
        Next
    End Function
    <Extension()>
    Public Function GenBankMetaData(Seq As Bio.ISequence) As Bio.IO.GenBank.GenBankMetadata
        If Seq.Metadata.ContainsKey(Bio.Util.Helper.GenBankMetadataKey) Then
            Return (Seq.Metadata(Bio.Util.Helper.GenBankMetadataKey))
        Else
            Return Nothing
        End If
    End Function

    <Extension()>
    Public Function HasGenBankMetadata(Seq As Bio.ISequence) As Boolean
        If IsNothing(Seq.Metadata) = True Then
            Return False
        Else
            Return True
        End If
    End Function
    <Extension()>
    Public Iterator Function HasGenBankMetadata(seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of Boolean)
        For Each Seq In seqs
            Yield Seq.HasGenBankMetadata
        Next
    End Function
    <Extension()>
    Public Function HasFeatures(Seq As Bio.ISequence) As Boolean
        If IsNothing(Seq.GenBankMetaData) = True Then
            Return False
        Else
            If IsNothing(Seq.GenBankMetaData.Features) = True Then
                Return False
            Else
                Return True
            End If
        End If
    End Function
    <Extension()>
    Public Iterator Function HasFeatures(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of Boolean)
        For Each Seq In Seqs
            Yield Seq.HasFeatures
        Next
    End Function
#End Region

#Region "Locus"
    ''' <summary>
    ''' Return LocusName or Empty String
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function LocusName(Seq As Bio.ISequence) As String
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Locus) = True OrElse IsNothing(md.Locus.Name) = True Then Return Nothing
        Return md.Locus.Name
    End Function
    ''' <summary>
    ''' Return Locus Date or Nothing
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function [Date](Seq As Bio.ISequence) As Date
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Locus) = True OrElse IsNothing(md.Locus.Date) = True Then Return Nothing
        Return md.Locus.Date
    End Function
    ''' <summary>
    ''' Return Locus Sequence Division Code or Nothing
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function SequenceDivisionCode(Seq As Bio.ISequence) As Bio.IO.GenBank.SequenceDivisionCode
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Locus) = True OrElse IsNothing(md.Locus.DivisionCode) = True Then Return Nothing
        Return md.Locus.DivisionCode
    End Function
    ''' <summary>
    ''' Return Locus MoleculeType Code or Nothing
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function MoleculeType(Seq As Bio.ISequence) As Bio.IO.GenBank.MoleculeType
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Locus) = True OrElse IsNothing(md.Locus.MoleculeType) = True Then Return Nothing
        Return md.Locus.MoleculeType
    End Function
    ''' <summary>
    ''' Return Locus SequenceType  or empty String
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function SequenceType(Seq As Bio.ISequence) As String
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Locus) = True OrElse IsNothing(md.Locus.SequenceType) = True Then Return Nothing
        Return md.Locus.SequenceType
    End Function
    ''' <summary>
    ''' Return Locus SequenceStrandType  or Nothing
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function SequenceStrandType(Seq As Bio.ISequence) As Bio.IO.GenBank.SequenceStrandType
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Locus) = True OrElse IsNothing(md.Locus.Strand) = True Then Return Nothing
        Return md.Locus.Strand
    End Function
    ''' <summary>
    ''' Return Locus SequenceStrandTopology  or Nothing
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function SequenceStrandTopology(Seq As Bio.ISequence) As Bio.IO.GenBank.SequenceStrandTopology
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True OrElse IsNothing(md.Locus) = True OrElse IsNothing(md.Locus.Strand) = True Then Return Nothing
        Return md.Locus.StrandTopology
    End Function
    ''' <summary>
    ''' Return LocusInfo or Nothing
    ''' </summary>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function Locus(Seq As Bio.ISequence) As Bio.IO.GenBank.GenBankLocusInfo
        Dim md = Seq.GenBankMetaData
        If IsNothing(md) = True Then Return Nothing
        Return md.Locus
    End Function

#End Region

#Region "Translate"

    <Extension()>
    Public Function Translate(Seq As Bio.ISequence, Optional Frame As Integer = 1, Optional RemoveTermianlStopCOdon As Boolean = True) As Bio.ISequence

        Dim RNASeq As Bio.ISequence
        If Frame > 0 Then
            Frame -= 1
            RNASeq = Bio.Algorithms.Translation.Transcription.Transcribe(Seq)

        Else
            Frame = System.Math.Abs(Frame + 1)
            RNASeq = Bio.Algorithms.Translation.Transcription.Transcribe(Seq.GetReverseComplementedSequence)
        End If
        Dim AASeq = Bio.Algorithms.Translation.ProteinTranslation.Translate(RNASeq, Frame)
        AASeq.ID = Seq.ID
        If RemoveTermianlStopCOdon = True Then
            Return AASeq.RemoveTerminalStop
        Else
            Return AASeq
        End If
    End Function
    <Extension()>
    Public Iterator Function Translate(Seqs As IEnumerable(Of Bio.ISequence), Optional Frame As Integer = 1, Optional RemoveTermianlStopCodon As Boolean = True) As IEnumerable(Of Bio.ISequence)

        For Each Seq In Seqs
            Yield Seq.Translate(Frame, RemoveTermianlStopCodon)
        Next
    End Function
#End Region

#Region "Unique Duplicates"

    <Extension()>
    Public Iterator Function Firsts(Seqs As IEnumerable(Of List(Of Bio.ISequence))) As IEnumerable(Of Bio.ISequence)
        For Each Seq In Seqs
            Yield Seq.First
        Next
    End Function

    <Extension()>
    Public Iterator Function Distinct_ByID(SeqsA As List(Of Bio.ISequence), SeqsB As List(Of Bio.ISequence)) As IEnumerable(Of Bio.ISequence)
        Dim c As New Szunyi.Sequences.Sorters.ByID
        For Each Seq In SeqsA
            If SeqsB.BinarySearch(Seq, c) < 0 Then
                Yield Seq
            End If
        Next
    End Function
    <Extension()>
    Public Iterator Function Intersect_ByID(SeqsA As List(Of Bio.ISequence), SeqsB As List(Of Bio.ISequence)) As IEnumerable(Of Bio.ISequence)
        Dim c As New Szunyi.Sequences.Sorters.ByID
        For Each Seq In SeqsA
            If SeqsB.BinarySearch(Seq, c) > -1 Then
                Yield Seq
            End If
        Next
    End Function
#End Region

#Region "Rename"
    <Extension()>
    Public Sub ReName_Increment_at_End(Seqs As IEnumerable(Of Bio.ISequence), Optional Separator As String = "_")
        Dim Index As Integer = 0

        For Each Seq In Seqs
            Seq.ID = Seq.ID & Separator & Index
            Index += 1
        Next

    End Sub
    <Extension()>
    Public Sub ReName_Increment_at_Start(Seqs As List(Of Bio.ISequence), Optional Separator As String = "_")
        Dim Index As Integer = 0
        For Each Seq In Seqs
            Seq.ID = Index & Separator & Seq.ID
            Index += 1
        Next

    End Sub
#End Region

#Region "BinarySearch"
    ''' <summary>
    ''' Return Seq with Same ID or Nothing 
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <param name="Seq"></param>
    ''' <returns></returns>
    <Extension()>
    Public Function SearchByID(Seqs As List(Of Bio.ISequence), Seq As Bio.ISequence) As Bio.ISequence
        Dim Index = Seqs.BinarySearch(Seq, SequnceID_Sorter)
        If Index > -1 Then
            Return Seqs(Index)
        Else
            Return Nothing
        End If
    End Function
#End Region

#Region "Clone, Convert"
    <Extension()>
    Public Function Clone(Seq As Bio.ISequence)
        If IsNothing(Seq) = True Then Return Nothing
        If Seq.GetType.Name = "QualitativeSequence" Then
            Dim s As QualitativeSequence = Seq
            Dim x As New Bio.QualitativeSequence(Seq.Alphabet, s.FormatType, Seq.ToArray, s.GetEncodedQualityScores)
            x.ID = Seq.ID
            Return x
        Else
            Dim x As New Bio.Sequence(Seq.Alphabet, Seq.ToArray)
            x.ID = Seq.ID
            If IsNothing(Seq.GenBankMetaData) = False Then
                x.Metadata.Add(Bio.Util.Helper.GenBankMetadataKey, Seq.GenBankMetaData.Clone)
            End If
            Return x
        End If

    End Function
    <Extension()>
    Public Iterator Function Clone(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of Bio.ISequence)
        For Each Seq In Seqs
            Yield Seq.Clone
        Next
    End Function
    <Extension()>
    Public Function ConvertTo_BioSeq(Seq As Bio.ISequence) As Bio.Sequence
        If Seq.GetType().Name = "QualitativeSequence" Then
            Dim t1 As Bio.QualitativeSequence = Seq
            Dim t2 As New Bio.Sequence(Alphabets.AmbiguousDNA, t1.ToArray)
            t2.ID = t1.ID
            Return t2
        Else
            Return Seq
        End If
    End Function
    <Extension()>
    Public Iterator Function ConvertTo_BioSeq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of Bio.Sequence)
        For Each Seq In Seqs
            Yield Seq.ConvertTo_BioSeq
        Next
    End Function
    <Extension()>
    Public Sub ConvertTo_GenBank(ByRef seqs As IEnumerable(Of ISequence))
        For Each Seq In seqs
            If Seq.Metadata.ContainsKey(Bio.Util.Helper.GenBankMetadataKey) = False Then
                Seq.GenBankMetaData_NACreate
            End If
        Next
    End Sub
#End Region

#Region "Group"

    <Extension()>
    Public Iterator Function GroupBy_ID(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))
        Dim x = From t In Seqs Select t Group By t.ID Into Group

        For Each g In x
            Yield g.Group.ToList
        Next
    End Function
    <Extension()> Public Iterator Function GroupBy_Seq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))
        Dim x = From t In Seqs Select t Group By t.ToArray Into Group

        For Each g In x
            Yield g.Group.ToList
        Next
    End Function
    <Extension()>
    Public Iterator Function GroupBy_ID_Seq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))
        Dim x = From t In Seqs Select t Group By t.ToArray, t.ID Into Group

        For Each g In x
            Yield g.Group.ToList
        Next
    End Function
#Region "Duplicates"
    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function DuplicatesBy_ID(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))

        For Each g In Seqs.GroupBy_ID
            If g.Count > 1 Then
                Yield g
            End If
        Next

    End Function

    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function DuplicatesBy_Seq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))

        For Each g In Seqs.GroupBy_Seq
            If g.Count > 1 Then
                Yield g
            End If
        Next

    End Function

    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function DuplicatesBy_ID_Seq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))

        For Each g In Seqs.GroupBy_ID_Seq
            If g.Count > 1 Then
                Yield g
            End If
        Next

    End Function



#End Region

#Region "Unique"
    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function UniqueBy_ID(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))

        For Each g In Seqs.GroupBy_ID
            If g.Count = 1 Then
                Yield g
            End If
        Next

    End Function

    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function UniqueBy_Seq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))

        For Each g In Seqs.GroupBy_Seq
            If g.Count = 1 Then
                Yield g
            End If
        Next

    End Function

    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function UniqueBy_ID_Seq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of List(Of Bio.ISequence))

        For Each g In Seqs.GroupBy_ID_Seq
            If g.Count = 1 Then
                Yield g
            End If
        Next

    End Function
#End Region

#Region "OneCopy"
    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function OneCopyBy_ID(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of Bio.ISequence)

        For Each g In Seqs.GroupBy_ID
            Yield g.First
        Next

    End Function

    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function OneCopyBy_Seq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of Bio.ISequence)

        For Each g In Seqs.GroupBy_Seq
            Yield g.First
        Next

    End Function

    ''' <summary>
    ''' It is Return All Of The Duplicated Sequence
    ''' </summary>
    ''' <param name="Seqs"></param>
    ''' <returns></returns>
    <Extension()>
    Public Iterator Function OneCopyBy_ID_Seq(Seqs As IEnumerable(Of Bio.ISequence)) As IEnumerable(Of Bio.ISequence)

        For Each g In Seqs.GroupBy_ID_Seq
            Yield g.First
        Next

    End Function

#End Region
#End Region


End Module