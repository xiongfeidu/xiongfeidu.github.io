#!/bin/bash

# sleep 9h

while [ true ];
do
	cd ~/Desktop/Elections\ Analysis/Java/src/main/resources/web
	python web2016.py 2026governors
	# python web2016.py 2021primaries-governors
	python web2016.py 2026house
	python web2016.py 2026senate
	# python web2016.py 2024president
	# python web2016.py 2024repprimaries
	# python web2016.py 2024demprimaries
	cd ~/xiongfeidu.github.io/election/results
	cp -r ../../../Desktop/Elections\ Analysis/Java/src/main/resources/web/2026 .
	# cp -r ../../../Desktop/Elections\ Analysis/Java/src/main/resources/web/2025/governors 2025
	# cp -r ../../../Desktop/Elections\ Analysis/Java/src/main/resources/web/2024 .
	DATE="$(date +"%F")"
	TIME="$(date +"%R")"
	#echo "update ${DATE} - update 2020 election results ${TIME}"
	git add 2026
	git commit -m "update ${DATE} - update 2026 election results ${TIME} EST"
	git push
	# sleep 59m
	# sleep 51s
	sleep 5m
	sleep 51s
done