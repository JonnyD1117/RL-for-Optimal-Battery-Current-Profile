from gym.envs.registration import register

register(
    id='battery-v0',
    entry_point='gym_battery.envs:BatteryEnv',
)