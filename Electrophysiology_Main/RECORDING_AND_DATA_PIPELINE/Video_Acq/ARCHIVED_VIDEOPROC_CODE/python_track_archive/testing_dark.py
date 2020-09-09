# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
ii=0
thiscont = contours[ii]
# find area, perimeter and circularity of each contour
area = cv2.contourArea(thiscont)
perimeter = cv2.arcLength(thiscont,True)
circ = 4*math.pi*(area/(perimeter**2))
# compute moments and find center of each contour
Mii = cv2.moments(thiscont)
x = int(Mii['m10']/Mii['m00'])
y = int(Mii['m01']/Mii['m00'])
center[0].append(x)
center[1].append(y)
# find rotated bounding rectangle for blob and get dimensions
minrect = cv2.minAreaRect(thiscont)
w_rect = minrect[1][0]
h_rect = minrect[1][1]
lin = w_rect/h_rect
# check linearity threshold
if lin > linThresh or lin < 1./linThresh:
    tooLinear = True
#                    print 'too linear!'
else:
    tooLinear = False
#                tCorner = tuple([int(i) for i in list(minrect[0])])
#                bCorner = tuple([int(i) for i in list(tuple(np.subtract(topCorner,minrect[1])))])

# find bounding rectangle to check area fraction
_,_,w_box,h_box = cv2.boundingRect(thiscont)
box_area = w_box*h_box
area_fraction = area/box_area

# find distance from largest blob
# obviously, this will be 0 for ii = 0
largest_center = np.array([center[0][0],center[1][0]])
this_center = np.array([center[0][ii],center[1][ii]])
dist_from_largest = np.linalg.norm(largest_center-this_center)

print 'Area: {}'.format(area)
print 'Area fraction: {}'.format(area_fraction)
print 'Box area: {}'.format(box_area)
print 'Too linear?: {}'.format(tooLinear)
print 'W_rect: {} ; h_rect: {}. Linearity: {}'.format(w_rect,h_rect,lin)
print 'Circularity: {}'.format(circ)
print 'Far out-ness: {}'.format(dist_from_largest)

cv2.destroyAllWindows()


testframe = colorframe.copy()


cv2.drawContours(testframe,[contours[ii]],0,(0,0,255),-1)
cv2.circle(testframe, (mcx,mcy), 5, (255,0,0), -1)

cv2.imshow('1',testframe)

