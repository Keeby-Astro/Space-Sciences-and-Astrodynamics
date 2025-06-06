% ALGORITHM 8.2: CALCULATION OF THE SPACECRAFT TRAJECTORY
% FROM PLANET 1 TO PLANET 2

function ...
[planet1, planet2, trajectory] = interplanetary(depart, arrive)
%{
    This function determines the spacecraft trajectory from the sphere
    of influence of planet 1 to that of planet 2 using Algorithm 8.2

    mu          - gravitational parameter of the sun (km^3/s^2)
    dum         - a dummy vector not required in this procedure

    planet_id   - planet identifier:
                  1 = Mercury
                  2 = Venus
                  3 = Earth
                  4 = Mars
                  5 = Jupiter
                  6 = Saturn
                  7 = Uranus
                  8 = Neptune
                  9 = Pluto

    year        - range: 1901 - 2099
    month       - range: 1 - 12
    day         - range: 1 - 31
    hour        - range: 0 - 23
    minute      - range: 0 - 60
    second      - range: 0 - 60

    jd1, jd2    - Julian day numbers at departure and arrival
    tof         - time of flight from planet 1 to planet 2 (s)
    Rp1, Vp1    - state vector of planet 1 at departure (km, km/s)
    Rp2, Vp2    - state vector of planet 2 at arrival (km, km/s)
    R1, V1      - heliocentric state vector of spacecraft at
                  departure (km, km/s)
    R2, V2      - heliocentric state vector of spacecraft at
                  arrival (km, km/s)

    depart      - [planet_id, year, month, day, hour, minute, second]
                  at departure
    arrive      - [planet_id, year, month, day, hour, minute, second]
                  at arrival

    planet1     - [Rp1, Vp1, jd1]
    planet2     - [Rp2, Vp2, jd2]
    trajectory  - [V1, V2]

    User M-functions required: planet_elements_and_sv, lambert
%}
% --------------------------------------------------------------------

global mu

planet_id = depart(1);
year      = depart(2);
month     = depart(3);
day       = depart(4);
hour      = depart(5);
minute    = depart(6);
second    = depart(7);

%...Use Algorithm 8.1 to obtain planet 1's state vector (don't
%...need its orbital elements ["dum"]):
[dum, Rp1, Vp1, jd1] = planet_elements_and_sv ...
              (planet_id, year, month, day, hour, minute, second);

planet_id = arrive(1);
year      = arrive(2);
month     = arrive(3);
day       = arrive(4);
hour      = arrive(5);
minute    = arrive(6);
second    = arrive(7);

%...Likewise use Algorithm 8.1 to obtain planet 2's state vector:
[dum, Rp2, Vp2, jd2] = planet_elements_and_sv ...
              (planet_id, year, month, day, hour, minute, second);

tof = (jd2 - jd1)*24*3600;

%...Patched conic assumption:
R1 = Rp1;
R2 = Rp2;

%...Use Algorithm 5.2 to find the spacecraft's velocity at
% departure and arrival, assuming a prograde trajectory:
[V1, V2] = lambert(R1, R2, tof, 'pro');

planet1    = [Rp1, Vp1, jd1];
planet2    = [Rp2, Vp2, jd2];
trajectory = [V1, V2];

end