// unique-distances.cc
// Solve puzzle to arrange markers on a grid such that all pairwise
// distances are unique.

// Written in 2020 by Scott McPeak, scott.g.mcpeak at gmail.com.
//
// To the extent possible under law, the author(s) have dedicated all
// copyright and related and neighboring rights to this software to the
// public domain worldwide. This software is distributed without any
// warranty.
//
// You should have received a copy of the CC0 Public Domain Dedication
// along with this software. If not, see
// <http://creativecommons.org/publicdomain/zero/1.0/>.

#include <algorithm>                   // sort, swap
#include <array>                       // array
#include <assert.h>                    // assert
#include <bitset>                      // bitset
#include <iostream>                    // cout, etc.
#include <list>                        // list
#include <stdint.h>                    // uint8_t

using namespace std;


// Return the maximum possible distance (squared) between two squares
// on the NxN checkerboard.
constexpr int MAX_DISTANCE(int N)
{
  return (N-1) * (N-1) * 2;
}


// An NxN checkerboard with up to N markers placed onto squares.
//
// This also keeps track of which of the possible distances have been
// used.
//
// Within this class, squares are numbered starting with 0 in the lower
// left, increasing by 1 as we move to the right, and increasing by N as
// we move up.  For example, the 3x3 board is numbered like this:
//
//   6 7 8
//   3 4 5
//   0 1 2
//
// Additionally, within this class, "distance" means the square of the
// cartesian distance, i.e., dx*dx + dy*dy.  The original statement of
// the puzzle asks for unique cartesian distances, but it is of course
// equivalent to instead use unique squares of distance, which avoids
// needlessly computing square roots.
//
template <int N>
class MarkedBoard {
public:      // class data
  // The tables here exist as a minor optimization for some commonly
  // executed arithmetic operations.  They are public because they are
  // not really specific to this class, but need to be class members
  // because it's not possible to have templatized global data.

  // Map from square number to its row and column, where the row is the
  // high 4 bits and the column is the low 4 bits.
  static uint8_t decomposeTable[N*N];

  // Map from 'sq1*N*N+sq2' to the distance between the squares.
  static uint16_t distanceTable[N*N*N*N];

private:     // instance data
  // Number of markers placed onto the board.
  //
  // Invariant: 0 <= m_numPlacedMarkers <= N
  //
  int m_numPlacedMarkers;

  // Set of square numbers that contain markers.  The meaningful entries
  // are those with indices in [0,m_numPlacedMarkers-1].  This array
  // records the order in which markers are placed, which is strictly
  // increasing to factor out the symmetry of marker placement order.
  //
  // Invariant: All meaningful are strictly increasing values in [0,N*N-1].
  //
  array<int,N> m_markers;

  // Set of distances between squares that have been used by some pair
  // of markers.  Not all distances are possible (only sums of squares
  // are), but I do not bother to compress away the unusable bits.
  //
  // Invariant: m_usedDistances.count() ==
  //            m_numPlacedMarkers * (m_numPlacedMarkers-1) / 2
  //
  bitset<MAX_DISTANCE(N) + 1> m_usedDistances;

public:      // methods
  // Construct an empty board.
  MarkedBoard()
    : m_numPlacedMarkers(0),
      m_markers(),                     // set below
      m_usedDistances()                // initially empty
  {
    // Initialize 'm_markers' to bogus values for determinism.
    for (int i=0; i < N; i++) {
      m_markers.at(i) = N;
    }
  }

  // Make a copy of the given board.
  MarkedBoard(MarkedBoard const &obj)
    : m_numPlacedMarkers(obj.m_numPlacedMarkers),
      m_markers         (obj.m_markers),
      m_usedDistances   (obj.m_usedDistances)
  {}

  // Return number of markers placed on board.
  //
  // Ensures: 0 <= return <= N
  //
  int numPlacedMarkers() const {
    return m_numPlacedMarkers;
  }

  // Return first square number where a marker could be added and
  // maintain the property that m_markers is in increasing order.  If
  // none is possible, this returns N*N.
  //
  // Ensures: 0 <= return <= N*N
  //
  int firstCandidateSquare() const {
    if (m_numPlacedMarkers == 0) {
      return 0;
    }
    else {
      // One greater than the last placed marker.
      return m_markers.at(m_numPlacedMarkers-1) + 1;
    }
  }

  // Return the square number for the marker at the given index.
  //
  // Requires: 0 <= markerIndex < m_numPlacedMarkers
  // Ensures: 0 <= return < N*N
  //
  int getSquareOfMarker(int markerIndex) const {
    return m_markers.at(markerIndex);
  }

  // True if the indicated square has a marker on it.
  //
  // Requires: 0 <= square < N*N
  //
  bool hasMarkerAt(int square) const {
    // This method is only used during printing, so it can be slow.
    for (int i=0; i < m_numPlacedMarkers; i++) {
      if (this->getSquareOfMarker(i) == square) {
        return true;
      }
    }
    return false;
  }

  // True if it appears we can place a marker at 'square'.  This might
  // return true in some cases where in fact we cannot place the marker;
  // that will be diagnosed by 'placeAt'.
  //
  // Requires: 0 <= square < N*N
  //
  bool canPlaceAt(int square) const {
    // This mechanism is no longer being used.
    return true;
  }

  // Modify this board, putting a marker on 'square'.  Return true iff
  // the marker can be placed without violating any constraints.  If
  // this returns false, the object is left in an indeterminate state
  // and should not be used.
  //
  // Requires: m_numPlacedMarkers < N
  // Requires: 0 <= square < N*N
  //
  bool placeAt(int square);

  // Return true if there is a pair of markers separated by the
  // given distance.
  //
  // Requires: 0 <= distance <= MAX_DISTANCE(N)
  //
  bool usedDistance(int distance) const {
    return m_usedDistances.test(distance);
  }

  // Return true if 'this' board is the same as 'other' except for some
  // possible rotation and/or reflection.
  bool isEquivalentTo(MarkedBoard const &other) const;
};


// Map from square number to its row and column, where the row is the
// high 4 bits and the column is the low 4 bits.
template <int N>
/*static*/ uint8_t MarkedBoard<N>::decomposeTable[N*N];

// Build 'decomposeTable'.
template <int N>
void makeDecomposeTable()
{
  assert(N <= 16);
  for (int r=0; r < N; r++) {
    for (int c=0; c < N; c++) {
      MarkedBoard<N>::decomposeTable[r*N + c] = (r << 4) + c;
    }
  }
}


// Given a square number on an NxN board, return its row and column,
// numbered starting at 0 from the left and bottom.
//
// Requires: 0 <= sq < N*N
// Ensures: r == sq / N
// Ensures: c == sq % N
//
template <int N>
void decompose(int &r, int &c, int sq)
{
  assert(0 <= sq && sq < N*N);

  // Use lookup table for around 10% better performance (when N==10 at
  // least) than the obvious implementation.
  uint8_t rc = MarkedBoard<N>::decomposeTable[sq];
  r = rc >> 4;
  c = rc & 0xf;

  // Check equivalence with slower method (if not NDEBUG).
  assert(r == sq / N);
  assert(c == sq % N);
}


// Given two square numbers on an NxN board, return the square of the
// cartesian distance between them.
//
// Requires: 0 <= sq1 < N*N
// Requires: 0 <= sq2 < N*N
// Ensures: 0 <= return <= (N-1)*(N-1)*2
//
template <int N>
int simpleComputeDistance(int sq1, int sq2)
{
  int r1, c1;
  decompose<N>(r1, c1, sq1);

  int r2, c2;
  decompose<N>(r2, c2, sq2);

  int dr = r1 - r2;
  int dc = c1 - c2;
  return dr*dr + dc*dc;
}


template <int N>
/*static*/ uint16_t MarkedBoard<N>::distanceTable[N*N*N*N];

template <int N>
void makeDistanceTable()
{
  for (int sq1=0; sq1 < N*N; sq1++) {
    for (int sq2=0; sq2 < N*N; sq2++) {
      int d = simpleComputeDistance<N>(sq1, sq2);
      MarkedBoard<N>::distanceTable[sq1 * N*N + sq2] = d;

      // Check that nothing was truncated.
      assert(MarkedBoard<N>::distanceTable[sq1 * N*N + sq2] == d);
    }
  }
}


// Same specification as 'simpleComputeDistance', but about 10% faster.
template <int N>
int computeDistance(int sq1, int sq2)
{
  int ret = MarkedBoard<N>::distanceTable[sq1 * N*N + sq2];
  assert(ret == simpleComputeDistance<N>(sq1, sq2));
  return ret;
}


template <int N>
bool MarkedBoard<N>::placeAt(int square)
{
  // Cannot place more than N markers.
  if (m_numPlacedMarkers >= N) {
    return false;
  }

  // Must place markers in increasing order.
  if (square < firstCandidateSquare()) {
    return false;
  }

  // Add the distance from the new square to each of the existing
  // markers.
  for (int i=0; i < m_numPlacedMarkers; i++) {
    int existing = m_markers.at(i);
    int distance = computeDistance<N>(square, existing);
    if (m_usedDistances.test(distance)) {
      // Distance is not unique.
      return false;
    }
    m_usedDistances.set(distance);
  }

  // Place the new marker.
  m_markers.at(m_numPlacedMarkers++) = square;

  return true;
}


// Return the square number corresponding to 'square' after applying
// transformations 'flip' and 'rot'.
//
// Requires: 0 <= square < N*N
// Requires: 0 <= flip < 2
// Requires: 0 <= rot < 4
//
template <int N>
int transformSquare(int square, int flip, int rot)
{
  int r, c;
  decompose<N>(r, c, square);

  if (flip) {
    // Flip across the main diagonal.
    swap(r, c);
  }

  // Rotate 90 degrees counterclockwise 'rot' times.
  while (rot--) {
    // Distance from left (c) becomes new distance from bottom (nr).
    int nr = c;

    // Distance from bottom (r) becomes new distance from right ((N-1)-nc).
    int nc = (N-1)-r;

    r = nr;
    c = nc;
  }

  return r*N + c;
}


template <int N>
bool MarkedBoard<N>::isEquivalentTo(MarkedBoard<N> const &other) const
{
  if (this->numPlacedMarkers() != other.numPlacedMarkers()) {
    return false;
  }

  // They must have the same set of used distances.
  for (int distance = 0; distance <= MAX_DISTANCE(N); distance++) {
    if (this->usedDistance(distance) != other.usedDistance(distance)) {
      return false;
    }
  }

  // Consider all 8 possible transformations.
  for (int flip=0; flip<2; flip++) {
    for (int rot=0; rot<4; rot++) {
      // Get the squares of 'other' after applying the transformation.
      array<int, N> otherSquares;
      for (int i=0; i < other.numPlacedMarkers(); i++) {
        int sq = other.getSquareOfMarker(i);
        int sqTrans = transformSquare<N>(sq, flip, rot);
        otherSquares.at(i) = sqTrans;
      }

      // Sort them to enable direct comparison to 'this'.
      sort(begin(otherSquares), end(otherSquares));

      // Compare the sets of squares, taking advantage of the fact that
      // both are in increasing order.
      for (int i=0; i < this->numPlacedMarkers(); i++) {
        int thisSquare = this->getSquareOfMarker(i);
        int otherSquare = otherSquares.at(i);
        if (thisSquare != otherSquare) {
          goto squareIsDifferent;
        }
      }

      // All squares were the same after transformation, so the boards
      // are equivalent.
      return true;

    squareIsDifferent:
      // Go on to try remaining transformations.
      ;
    }
  }

  // No transformation made them equivalent.
  return false;
}


// Return true if 'square' is on the main diagonal of the NxN board.
template <int N>
bool onDiagonal(int square)
{
  int r, c;
  decompose<N>(r, c, square);

  return r == c;
}


// Return true if 'square' is strictly above the main diagonal of the
// NxN board.
template <int N>
bool aboveDiagonal(int square)
{
  int r, c;
  decompose<N>(r, c, square);

  return r > c;
}


// Return true if 'square' is in the south-southwest octant of the
// NxN board.  For example, in the 4x4 case, the 'X' letters indicate
// that octant:
//
//   r:
//   3  ....
//   2  ....
//   1  .X..
//   0  XX..
//   c: 0123
//
// and for the 3x3 case:
//
//   r:
//   2  ...
//   1  .X.
//   0  XX.
//   c: 012
//
template <int N>
bool inSSWOctant(int square)
{
  if (aboveDiagonal<N>(square)) {
    return false;
  }

  int r, c;
  decompose<N>(r, c, square);
  if (c >= (N+1)/2) {
    return false;            // Right of midline.
  }

  return true;
}


// Return true if all placed markers on 'board' are on the main diagonal.
template <int N>
bool allMarkersOnDiagonal(MarkedBoard<N> const &board)
{
  for (int i=0; i < board.numPlacedMarkers(); i++) {
    if (!onDiagonal<N>(board.getSquareOfMarker(i))) {
      return false;          // Marker is not on diagonal.
    }
  }
  return true;
}


// Number of calls to 'placeMarkers', kept in order to measure the
// effectiveness of some optimizations.
static long numPlaceMarkersCalls = 0;


// Add to 'solutions' all of the solved boards obtainable by adding
// markers to 'orig' and that are not already present, including as
// rotated and/or flipped variants, in 'solutions'.
//
// If 'orig' is already solved, print it alone.
//
template <int N>
void placeMarkers(list<MarkedBoard<N> > &solutions,
                  MarkedBoard<N> const &orig)
{
  numPlaceMarkersCalls++;

  if (orig.numPlacedMarkers() == N) {
    // Found a solution.  Check it for uniqueness.
    for (auto const &sol : solutions) {
      if (orig.isEquivalentTo(sol)) {
        return;
      }
    }

    solutions.push_back(orig);
    return;
  }

  for (int square = orig.firstCandidateSquare(); square < N*N; square++) {
    // Require that the first marker be placed in the SSW octant.  That
    // reduces the number of symmetric configurations that are produced
    // by the enumeration before explicitly checking for symmetry.
    if (orig.numPlacedMarkers() == 0) {
      if (!inSSWOctant<N>(square)) {
        continue;
      }
    }

    // If all currently placed markers are on the diagonal, then only
    // place new markers on or below that diagonal, to further reduce
    // enumeration of symmetric configurations.
    else if (allMarkersOnDiagonal(orig) &&
             aboveDiagonal<N>(square)) {
      continue;
    }

    if (orig.canPlaceAt(square)) {
      MarkedBoard<N> modified(orig);
      if (modified.placeAt(square)) {
        placeMarkers(solutions, modified);
      }
    }
  }
}


// Print a 2D diagram of the board with markers, like:
//
// Label:
//   OO.
//   ..O
//   ...
//
// for one of the 3x3 solutions.
template <int N>
void printBoard(ostream &os, string const &label, MarkedBoard<N> const &board)
{
  os << label << '\n';

  for (int row = N-1; row >= 0; row--) {
    os << "  ";
    for (int col = 0; col < N; col++) {
      int square = row*N + col;
      os << (board.hasMarkerAt(square)? 'O' : '.');
    }
    os << '\n';
  }

  os << "  Marked squares:";
  for (int i=0; i < board.numPlacedMarkers(); i++) {
    os << ' ' << board.getSquareOfMarker(i);
  }
  os << '\n';

  os << "  Used distances:";
  for (int distance = 0; distance <= MAX_DISTANCE(N); distance++) {
    if (board.usedDistance(distance)) {
      os << ' ' << distance;
    }
  }
  os << '\n';
}


// Return the number of unique distances between pairs of distinct
// squares on the NxN board.
template <int N>
int numUniqueDistances()
{
  bitset<MAX_DISTANCE(N)+1> distances;

  for (int a = 0; a < N*N - 1; a++) {
    if (aboveDiagonal<N>(a)) {
      continue;
    }
    for (int b = a+1; b < N*N; b++) {
      if (aboveDiagonal<N>(b)) {
        continue;
      }
      distances.set(computeDistance<N>(a, b));
    }
  }

  return distances.count();
}


// Print all unique (up to rotation and reflection) solutions for the
// NxN board.
template <int N>
void printSolutions()
{
  makeDecomposeTable<N>();
  makeDistanceTable<N>();

  int numPairs = (N * (N-1) / 2);
  int numDistances = numUniqueDistances<N>();

  cout << "N=" << N
       << " pairs=" << numPairs
       << " distances=" << numDistances
       << '\n';

  if (numPairs > numDistances) {
    // At N=16 and beyond, the number of distinct pairs is more than the
    // number of possible distances, so it's clearly impossible.
    //
    // This happens because each Pythagorean triple reuses a distance
    // that can be achieved along a single row or column, so as the
    // number of triples increases, eventually there aren't enough
    // unique distances to place all markers.
    cout << "No solutions exist because numPairs exceeds numDistances.\n";
    return;
  }

  list<MarkedBoard<N> > solutions;
  MarkedBoard<N> emptyBoard;
  placeMarkers(solutions, emptyBoard);

  cout << "numPlaceMarkersCalls: " << numPlaceMarkersCalls << '\n';
  cout << "Found " << solutions.size() << " solutions.\n";
  for (auto const &sol : solutions) {
    printBoard(cout, "Solution:", sol);
  }
}


int main()
{
  printSolutions<10>();
  return 0;
}


// EOF
