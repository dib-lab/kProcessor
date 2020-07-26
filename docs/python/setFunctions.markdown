# kDataFrame Set Functions

Set functions are functions that performs Union, Intersection and Difference between kDataFrames.

> **Supports only kDataFrameMQF**

## Parameters

All set functions takes list of kDataFrames as input and returns a new kDataFrame after performing the set function


## Union

```python

union_kf = kp.kFrameUnion([kDataFrame1, kDataFrame2])

```
## Intersect

```python

intersect_kf = kp.kFrameIntersect([kDataFrame1, kDataFrame2])

```
## Difference

```python

diff_kf = kp.kFrameDiff([kDataFrame1, kDataFrame2])

```
