function AnnounceStartVid(obj, event, string_arg)

event_time = datestr(event.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF');

msg = [string_arg ' at ' event_time];
disp(msg);
