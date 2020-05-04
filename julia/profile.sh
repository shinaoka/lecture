#!/bin/sh

echo begin >> test_profile
date >> test_profile

date > profile
julia profile_loop_update.jl profile.ini 1>> profile 2>> test_profile
date >> profile

echo end >> test_profile
date >> test_profile
