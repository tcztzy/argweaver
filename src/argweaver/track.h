#ifndef ARGWEAVER_TRACK_H
#define ARGWEAVER_TRACK_H

// c++ includes
#include <string.h>
#include <string>
#include <vector>

// arghmm includes
#include "logging.h"
#include "parsing.h"


namespace argweaver {

using namespace std;

// A region within a chromosome
// start is inclusive, end is exclusive
class Region
{
public:
    Region(string chrom="", int start=0, int end=0) :
        chrom(chrom), start(start), end(end)
    {}

    void set(const string &_chrom, int _start, int _end) {
        chrom = _chrom;
        start = _start;
        end = _end;
    }

    int length() const {
        return end - start;
    }

    string chrom;
    int start;
    int end;
};


// A region within a chromosome associated with a value
// start is inclusive, end is exclusive
template <class T>
class RegionValue {
public:
    RegionValue() :
        chrom(""), start(0), end(0)
    {}

    RegionValue(string chrom, int start, int end, T value) :
        chrom(chrom), start(start), end(end), value(value)
    {}

    int length() const {
        return end - start;
    }

    string chrom;
    int start;
    int end;
    T value;
};


typedef char NullValue;
typedef RegionValue<char> RegionNullValue;


// A track is a series of regions each associated with a value
template <class T>
class Track : public vector<RegionValue<T> > {
protected:
    //used for the find function
    int start_idx;
public:


    // Returns start coordinate if regions are available
    // Returns -1 otherwise
    int start_coord() const {
        if (Track<T>::size() == 0)
            return -1;
        else
            return Track<T>::front().start;
    }

    // Returns end coordinate if regions are available
    // Returns -1 otherwise
    int end_coord() const {
        if (Track<T>::size() == 0)
            return -1;
        else
            return Track<T>::back().end;
    }

    // Returns value of region containing position
    // If start_idx not NULL, is updated to index of return value
    T find(int pos, const T &default_value, int *start_idx=NULL) const {
        int start = ( start_idx == NULL ? 0 : *start_idx );
        for (unsigned int i=start; i<Track<T>::size(); i++) {
            const RegionValue<T> &region = Track<T>::at(i);
            if (region.start <= pos && pos < region.end) {
                if (start_idx != NULL) *start_idx = i;
                return region.value;
            }
        }
        // region not found
        if (start_idx != NULL) *start_idx = 0;
        return default_value;
    }

    // Returns index of region containing position pos
    int index(int pos) const {
        for (unsigned int i=0; i<Track<T>::size(); i++) {
            const RegionValue<T> &region = Track<T>::at(i);
            if (region.start <= pos && pos < region.end)
                return i;
        }
        // region not found
        return -1;
    }

    // Adds one region to the track
    void append(string chrom, int start, int end, T value) {
        this->push_back(RegionValue<T>(chrom, start, end, value));
    }

    // Reads one region from a map file and adds it to the track
    bool read_track_line(const char *line)
    {
        string chrom;
        int start, end;
        T value;

        if (!read_track_line(line, chrom, start, end, value))
            return false;
        append(chrom, start, end, value);
        return true;
    }

    // combines adjacent entries in sorted track
    // loses "value" term; mainly intended for TrackNullValue
    void merge() {
        Track<T> oldvec = *this;
        this->clear();
        for (unsigned int i=0; i < oldvec.size(); ) {
            unsigned int j=i+1;
            int currEnd = oldvec[i].end;
            const string chrom = oldvec[i].chrom;
            while (j < oldvec.size()) {
                if (oldvec[j].chrom == chrom &&
                    oldvec[j].start <= currEnd) {
                    currEnd = oldvec[j].end;
                    j++;
                } else break;
            }
            this->push_back(RegionValue<T>(oldvec[i].chrom,
                                           oldvec[i].start, currEnd,
                                           oldvec[i].value));
            i = j;
        }
    }


    // adds entries in track t to this track, assuming both are sorted
    // note: Does not properly merge values; meant for use with TrackNullValue
    // also: requires that all tracks come from a single chromosome
    void merge_tracks(const Track<T> t) {
        Track<T> oldvec = *this;
        this->clear();
        unsigned int idx0=0, idx1=0;
        if (t.size() == 0 && oldvec.size() == 0) return;
        const string chrom = ( t.size() == 0 ? oldvec[0].chrom : t[0].chrom );
        if (t.size() > 0)
        while (idx1 < t.size() && idx0 < oldvec.size()) {
            const RegionValue<T> &region0 = oldvec[idx0];
            const RegionValue<T> &region1 = t[idx1];
            assert(region0.chrom == chrom);
            assert(region1.chrom == chrom);
            if (region0.end < region1.start) {
                this->push_back(region0);
                idx0++;
            } else if (region1.end < region0.start) {
                this->push_back(region1);
                idx1++;
            } else {
                this->push_back(RegionValue<T>(region0.chrom,
                                               min(region0.start, region1.start),
                                               max(region0.end, region1.end),
                                               region0.value));
                idx0++;
                idx1++;
            }
        }
        if (idx1 < t.size()) {
            assert(idx0 == oldvec.size());
            for ( ; idx1 < t.size(); idx1++) {
                assert(t[idx1].chrom == chrom);
                this->push_back(t[idx1]);
            }
        }
        if (idx0 < oldvec.size()) {
            assert(idx1 == t.size());
            for ( ; idx0 < oldvec.size(); idx0++) {
                assert(t[idx0].chrom == chrom);
                this->push_back(oldvec[idx0]);
            }
        }
        this->merge();
    }
};


typedef Track<NullValue> TrackNullValue;

// Reads one region from a map file
template <class T>
bool read_track_line(const char *line, RegionValue<T> &region);


// A reader for reading a track from a map file
template <class T>
class TrackReader {
public:
    TrackReader() :
        has_error(false),
        line(NULL),
        linesize(1024),
        lineno(0)
    {}
    ~TrackReader() {
        if (line)
            delete [] line;
    }

    // open a map file by filename
    bool open(const char *filename) {
        FILE *infile;
        if ((infile = fopen(filename, "r")) == NULL) {
            printError("cannot read file '%s'", filename);
            return false;
        }

        lineno = 0;
        has_error = false;

        return true;
    }

    // open a map file by stream
    bool open(FILE *_infile) {
        infile = _infile;
        has_error = false;
        return true;
    }

    // Fetches the next RegionValue from a map file
    // Returns true if region read, false otherwise
    bool next(RegionValue<T> &region) {
        lineno++;
        if (fgetline(&line, &linesize, infile)) {
            // ignore track lines
            //if (!strncmp(line, "track", 5))
            //    continue;

            if (!read_track_line(line, region)) {
                // error reading line
                has_error = true;
                return false;
            }
        } else {
            // no more lines in file
            return false;
        }

        // line has been successfully read
        return true;
    }

    // Returns true if error encountered
    bool error() const {
        return has_error;
    }

    // Returns current line number in map file
    int line_number() const {
        return lineno;
    }


protected:
    bool has_error;
    FILE *infile;
    char *line;
    int linesize;
    int lineno;
};


// Reads a track from a map stream
template <class T>
bool read_track(FILE *infile, Track<T> *track)
{
    TrackReader<T> reader;
    RegionValue<T> region;
    reader.open(infile);

    while (reader.next(region)) {
        track->push_back(region);
    }
    if (reader.error()) {
        printError("could not read track line %d", reader.line_number());
        return false;
    }

    return true;
}


// Reads a track from a map stream
// NOTE: only regions within chrom:start-end are kept
template <class T>
bool read_track_filter(FILE *infile, Track<T> *track,
                       string chrom, int start, int end)
{
    TrackReader<T> reader;
    RegionValue<T> region;
    reader.open(infile);

    while (reader.next(region)) {
        // only keep regions that overlap desired region
        if (region.chrom == chrom &&
            region.end > start && region.start < end) {
            // trim region
            if (region.start < start)
                region.start = start;
            if (region.end > end)
                region.end = end;

            track->push_back(region);
        }
    }
    if (reader.error()) {
        printError("could not read track line %d", reader.line_number());
        return false;
    }

    return true;
}

template <class T>
bool read_track_filter(FILE *infile, Track<T> *track, const Region &region)
{
    return read_track_filter<T>(infile, track,
                                region.chrom, region.start, region.end);
}


// Reads a track from a map file
template <class T>
bool read_track(const char *filename, Track<T> *track)
{
    FILE *infile;
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return false;
    }

    bool result = read_track(infile, track);

    fclose(infile);
    return result;
}

// Reads a track from a map file
// NOTE: only regions with chrom:start-end are kept
template <class T>
bool read_track_filter(const char *filename, Track<T> *track,
                       string chrom, int start, int end)
{
    FILE *infile;
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return false;
    }

    bool result = read_track_filter(infile, track, chrom, start, end);

    fclose(infile);
    return result;
}



} // namespace argweaver

#endif // ARGWEAVER_TRACK
